import pandas as pd
import pysam
import re
import argparse
import os
# import gzip
import subprocess

def get_pairs(regions, chrom1, chrom2, break_pos1, break_pos2):
    """
    Fetch all paired reads, including those with split alignments

    Parameters:
    - regions: the two SV regions to collect reads from
    - chrom1: chromosome # of first region
    - chrom2: chromosome # of second region

    Returns:
    - split_alignments: contains all sets of paired reads with split alignments
    - paired_alignments: contains all sets of paired reads without split alignments
    """
    split_alignments = []
    nonsplit_alignments = []

    # Iterate over first breakpoint end
    for i, region in enumerate(regions):
        for aln in region:
            # If we've seen this read before, skip it
            if (
                aln.is_duplicate
                or any(aln in tup for tup in split_alignments)
                or any(aln in tup for tup in nonsplit_alignments)
            ):
                # print("Duplicate read or read already seen")
                continue

            if not (aln.is_mapped):
                # print("Unmapped read")
                continue

            try:
                mate = samfile.mate(aln)
            except ValueError as e:
                print(e)
                continue

            # Filter out unhelpful discordant pairs
            if not aln.is_proper_pair and not (
                aln.is_supplementary or aln.is_secondary or aln.has_tag("SA")
            ):
                if (
                    i == 0
                    and (
                        not mate.reference_name == chrom2
                        or not (
                            break_pos2 - 500 < mate.reference_start
                            and break_pos2 + 500 > mate.reference_start
                        )
                    )
                ) or (
                    i == 1
                    and (
                        not mate.reference_name == chrom1
                        or not (
                            break_pos1 - 500 < mate.reference_start
                            and break_pos1 + 500 > mate.reference_start
                        )
                    )
                ):
                    continue

            if aln.is_supplementary or aln.is_secondary or aln.has_tag("SA"):
                sa_tag = (aln.get_tag("SA")).split(";")[:-1]
                sa_tag = sa_tag[0]
                sup_name, sup_pos, sup_strand, sup_cigar, sup_mapq, sup_tlen = (
                    sa_tag.split(",")
                )
                sup_pos = int(sup_pos)
                sup_tlen = int(sup_tlen)

                found_split2 = False
                for sup_read in samfile.fetch(sup_name, sup_pos, sup_pos + 1):
                    if (
                        sup_read.query_name == aln.query_name
                        and sup_read.reference_name == sup_name
                        and (sup_read.reference_start + 1) == sup_pos
                        and re.sub(r"[SH]", "", sup_read.cigarstring)
                        == re.sub(r"[SH]", "", sup_cigar)
                    ):
                        non_sup = (
                            aln
                            if not aln.is_supplementary and not aln.is_secondary
                            else sup_read
                        )
                        sup = (
                            aln
                            if aln.is_supplementary or aln.is_secondary
                            else sup_read
                        )
                        read1 = non_sup if non_sup.is_read1 else mate
                        read2 = non_sup if non_sup.is_read2 else mate
                        split_alignments.append((read1, sup, read2))
                        found_split2 = True
                        break
                if not found_split2:
                    print(sup_name, sup_pos, sup_cigar, sup_tlen)
                    raise NameError("split mate not found")

            elif mate.is_supplementary or mate.is_secondary or mate.has_tag("SA"):
                sa_tag = (mate.get_tag("SA")).split(";")[:-1]
                sa_tag = sa_tag[0]
                sup_name, sup_pos, sup_strand, sup_cigar, sup_mapq, sup_tlen = (
                    sa_tag.split(",")
                )
                sup_pos = int(sup_pos)
                sup_tlen = int(sup_tlen)

                found_split2 = False
                for sup_read in samfile.fetch(sup_name, sup_pos, sup_pos + 1):
                    if (
                        sup_read.query_name == mate.query_name
                        and sup_read.reference_name == sup_name
                        and (sup_read.reference_start + 1) == sup_pos
                        and re.sub(r"[SH]", "", sup_read.cigarstring)
                        == re.sub(r"[SH]", "", sup_cigar)
                    ):
                        non_sup = (
                            mate
                            if not mate.is_supplementary and not mate.is_secondary
                            else sup_read
                        )
                        sup = (
                            mate
                            if mate.is_supplementary or mate.is_secondary
                            else sup_read
                        )
                        read1 = non_sup if non_sup.is_read1 else aln
                        read2 = non_sup if non_sup.is_read2 else aln
                        split_alignments.append((read1, sup, read2))
                        found_split2 = True
                        break
                if not found_split2:
                    print(sup_name, sup_pos, sup_cigar, sup_tlen)
                    raise NameError("split mate not found")

            elif not aln.is_proper_pair and aln.next_reference_name == (
                chrom2 if i == 0 else chrom1
            ):
                read1 = aln if aln.is_read1 else mate
                read2 = aln if aln.is_read2 else mate
                nonsplit_alignments.append((read1, read2))

            elif aln.is_proper_pair:
                read1 = aln if aln.is_read1 else mate
                read2 = aln if aln.is_read2 else mate
                nonsplit_alignments.append((read1, read2))

    return split_alignments, nonsplit_alignments


def fetch_alignments(
    bamfile,
    chrom1,
    pos1,
    chrom2,
    pos2,
    sv_type,
    read_support,
    features,
    orientation,
    hom_len,
    hom,
    region_size,
    mapq_threshold=15,
):
    """
    Fetch reads from a BAM file that have split alignments
    between two genomic coordinates as well as all paired-end reads.

    Parameters:
    - bamfile: str, path to the BAM file.
    - chrom1: str, chromosome for the first position.
    - pos1: int, first position in the genome.
    - chrom2: str, chromosome for the second position.
    - pos2: int, second position in the genome.
    - region_size: int, the size of the region around the positions to fetch alignments from.
    - mapq_threshold: int, cutoff for collecting reads above certain mapping quality.

    Returns:
    - split_df: DataFrame, details of split alignments.
    - paired_df: DataFrame, details of discordant alignments.
    """
    # Open the BAM file
    samfile = pysam.AlignmentFile(bamfile, "rb")

    # Fetch alignments region_size away from both positions
    region1 = samfile.fetch(chrom1, pos1 - region_size, pos1 + region_size + 1)
    region2 = samfile.fetch(chrom2, pos2 - region_size, pos2 + region_size + 1)

    split_alignments, nonsplit_alignments = get_pairs(
        [region1, region2], chrom1, chrom2, pos1, pos2
    )

    fast1_name = f"fastq/b_{chrom1}:{pos1+1}_{chrom2}:{pos2+1}_1.fastq"
    fast2_name = f"fastq/b_{chrom1}:{pos1+1}_{chrom2}:{pos2+1}_2.fastq"
    os.makedirs('fastq', exist_ok=True)
    with open(fast1_name, "w") as fast1:
        with open(fast2_name, "w") as fast2:
            for pair in nonsplit_alignments:
                read1, read2 = pair
                fast1.write(
                    read1.query_name
                    + "\n"
                    + read1.query_sequence
                    + "\n"
                    + "+\n"
                    + "".join(list(map(lambda x: chr(x + 33), read1.query_qualities)))
                    + "\n"
                )
                fast2.write(
                    read2.query_name
                    + "\n"
                    + read2.query_sequence
                    + "\n"
                    + "+\n"
                    + "".join(list(map(lambda x: chr(x + 33), read2.query_qualities)))
                    + "\n"
                )
            for pair in split_alignments:
                read1_1, read1_2, read2 = pair
                fast1.write(
                        read1_1.query_name
                        + "\n"
                        + read1_1.query_sequence
                        + "\n"
                        + "+\n"
                        + "".join(
                            list(map(lambda x: chr(x + 33), read1_1.query_qualities)))
                        + "\n"
                    )
                fast2.write(
                        read2.query_name
                        + "\n"
                        + read2.query_sequence
                        + "\n"
                        + "+\n"
                        + "".join(
                            list(map(lambda x: chr(x + 33), read2.query_qualities)))
                        + "\n"
                    )
    subprocess.run(['gzip', fast1_name])
    subprocess.run(['gzip', fast2_name])

    # Collect details of split alignments
    split_alignment_details = []
    for pair in split_alignments:
        read1_1, read1_2, read2 = pair

        # Filter out low quality reads
        if (
            (
                read1_1.mapping_quality > mapq_threshold
                and read1_2.mapping_quality > mapq_threshold
            )
            or read2.mapping_quality > mapq_threshold
        ) or ((read1_1.is_mapped and read1_2.is_mapped) or read2.is_mapped):
            for read in pair:
                split_alignment_details.append(
                    {
                        "break_chrom1": chrom1,
                        "break_pos1": pos1 + 1,
                        "break_chrom2": chrom2,
                        "break_pos2": pos2 + 1,
                        "break_sv_type": sv_type,
                        "break_read_support": read_support,
                        "break_features": features,
                        "break_orientation": orientation,
                        "AA_homology_len": hom_len,
                        "AA_homology_seq": hom,
                        "query_name": read.query_name,
                        "query_short": read.query_name.split(":")[-1],
                        "split": read.has_tag("SA"),
                        "proper_pair": (
                            "Concordant" if read.is_proper_pair else "Discordant"
                        ),
                        "read_num": "1" if read.is_read1 else "2",
                        "query_chrom": read.reference_name,
                        "query_pos": read.reference_start + 1,
                        "query_end": read.reference_end,
                        "query_orientation": "+" if read.is_forward else "-",
                        "query_cigar": read.cigarstring,
                        "query_aln_full": read.query_sequence,
                        "query_aln_sub": read.query_alignment_sequence,
                    }
                )

    # Collect details of nonsplit alignments
    nonsplit_alignment_details = []
    for pair in nonsplit_alignments:
        read1, read2 = pair
        # Filter out low quality reads
        if (
            read1.mapping_quality > mapq_threshold
            or read2.mapping_quality > mapq_threshold
        ) or (read1.is_mapped or read2.is_mapped):
            for read in pair:
                nonsplit_alignment_details.append(
                    {
                        "break_chrom1": chrom1,
                        "break_pos1": pos1 + 1,
                        "break_chrom2": chrom2,
                        "break_pos2": pos2 + 1,
                        "break_sv_type": sv_type,
                        "break_read_support": read_support,
                        "break_features": features,
                        "break_orientation": orientation,
                        "AA_homology_len": hom_len,
                        "AA_homology_seq": hom,
                        "query_name": read.query_name,
                        "query_short": read.query_name.split(":")[-1],
                        "split": read.has_tag("SA"),
                        "proper_pair": (
                            "Concordant" if read1.is_proper_pair else "Discordant"
                        ),
                        "read_num": "1" if read.is_read1 else "2",
                        "query_chrom": read.reference_name,
                        "query_pos": read.reference_start + 1,
                        "query_end": read.reference_end,
                        "query_orientation": "+" if read.is_forward else "-",
                        "query_cigar": read.cigarstring,
                        "query_aln_full": read.query_sequence,
                        "query_aln_sub": read.query_alignment_sequence,
                    }
                )

    samfile.close()

    # Convert details to pandas DataFrames
    split_df = pd.DataFrame(split_alignment_details)
    nonsplit_df = pd.DataFrame(nonsplit_alignment_details)
    # paired_df = pd.DataFrame(paired_alignment_details)

    return split_df, nonsplit_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collect reads from the ends of a SV for refinement."
    )
    parser.add_argument(
        "refine",
        help="Radius of refinement region",
        type=int,
        default=350,
    )
    parser.add_argument(
        "sum",
        help="path to AA amplicon summaries folder",
        type=str,
    )
    parser.add_argument(
        "bam",
        help="path to bamfile",
        type=str,
    )
    parser.add_argument(
        "-v", "--verbose", help="Include debug output", action="store_true"
    )
    parser.add_argument("-f", "--file", help="File to store alignments", type=str)
    args = parser.parse_args()
    # Define the BAM file and positions
    bamfile = args.bam
    samfile = pysam.AlignmentFile(bamfile, "rb")
    # print(samfile.header)
    df = pd.concat(
        [
            pd.read_csv(os.path.join(args.sum, f), sep="\t")
            for f in os.listdir(args.sum)
            if ".tsv" in f
        ],
        ignore_index=True,
    )
    # Define output table
    output = pd.DataFrame()
    num_split = 0
    num_discord = 0
    num_concord = 0
    num_paired = 0

    for index, row in df.iterrows():
        chrom1 = row["chrom1"]
        pos1 = row["pos1"]
        chrom2 = row["chrom2"]
        pos2 = row["pos2"]
        sv_type = row["sv_type"]
        read_support = row["read_support"]
        features = row["features"]
        orientation = row["orientation"]
        hom_len = row["homology_length"]
        homology = row["homology_sequence"]

        # Fetch the DataFrames for split alignments and paired alignments
        split_df, nonsplit_df = fetch_alignments(
            bamfile,
            chrom1,
            pos1,
            chrom2,
            pos2,
            sv_type,
            read_support,
            features,
            orientation,
            hom_len,
            homology,
            args.refine,
        )
        num_split += len(split_df) / 3
        num_paired += (2 * len(split_df) / 3) + len(nonsplit_df)
        # num_paired += len(paired_df)

        if args.verbose:
            print(f"Reads at junction defined by {chrom1} {pos1} and {chrom2} {pos2}:")
            print(f"Number split pairs: {len(split_df) / 3}")
            print(f"Number nonsplit pairs: {len(nonsplit_df) / 2}")
            print(f"Number paired: {(2 * len(split_df) / 3) + len(nonsplit_df)}")
        output = pd.concat([output, split_df], ignore_index=True)
        output = pd.concat([output, nonsplit_df], ignore_index=True)

    output.to_csv(
        args.file if args.file else args.bam.split("/")[-1].split(".")[0] + ".tsv",
        sep="\t",
        index=False,
    )

    if args.verbose:
        print("Total Number of split reads:", num_split)
        print("Total Number of paired reads:", num_paired)

    # df_nonsplit = pd.read_csv("alignments.tsv", sep="\t")
    # df_nonsplit = df_nonsplit[
    #     (df_nonsplit["split"] == False) & df_nonsplit["query_cigar"].str.contains("S")
    # ]
    # df_nonsplit.to_csv("special_cases.csv")
    # pd.set_option("display.max_colwidth", None)
    # df_nonsplit = df_nonsplit[
    #     [
    #         "split",
    #         "AA_homology_seq",
    #         "break_chrom1",
    #         "break_pos1",
    #         "break_chrom2",
    #         "break_pos2",
    #         "query_name",
    #         "query_chrom",
    #         "query_pos",
    #         "query_end",
    #         "query_orientation",
    #         "query_cigar",
    #         "query_aln_sub",
    #         "query_aln_full",
    #     ]
    # ]
    # df_nonsplit.query_name = df_nonsplit.query_name.apply(lambda x: x.split(":")[-1])
    # df_nonsplit
    # nonsplit_matches = df_nonsplit[
    #     (df_nonsplit["break_pos1"] == df_nonsplit["query_pos"])
    #     | (df_nonsplit["break_pos1"] == df_nonsplit["query_end"])
    #     | (df_nonsplit["break_pos2"] == df_nonsplit["query_pos"])
    #     | (df_nonsplit["break_pos2"] == df_nonsplit["query_end"])
    # ]

    # df_split = pd.read_csv("alignments.tsv", sep="\t")
    # df_split = df_split[df_split["split"] == True]
    # pd.set_option("display.max_colwidth", None)
    # df_split = df_split[
    #     [
    #         "split",
    #         "AA_homology_seq",
    #         "break_chrom1",
    #         "break_pos1",
    #         "break_chrom2",
    #         "break_pos2",
    #         "query_name",
    #         "query_chrom",
    #         "query_pos",
    #         "query_end",
    #         "query_orientation",
    #         "query_cigar",
    #         "query_aln_sub",
    #         "query_aln_full",
    #     ]
    # ]
    # df_split.query_name = df_split.query_name.apply(lambda x: x.split(":")[-1])
    # df_split
    # half_match = df_split[
    #     (df_split["break_pos1"] == df_split["query_pos"])
    #     | (df_split["break_pos1"] == df_split["query_end"])
    #     | (df_split["break_pos2"] == df_split["query_pos"])
    #     | (df_split["break_pos2"] == df_split["query_end"])
    # ]
    # full_match = half_match[
    #     half_match.groupby("query_name")["query_name"].transform("size") > 1
    # ]
    # full_match

    # new_awesomesauce_df = pd.concat([nonsplit_matches, full_match], axis=0).sort_values(
    #     ["AA_homology_seq", "split"]
    # )
    # new_awesomesauce_df["break_start"] = (
    #     new_awesomesauce_df["break_pos1"] == new_awesomesauce_df["query_pos"]
    # ) | (new_awesomesauce_df["break_pos2"] == new_awesomesauce_df["query_pos"])
    # new_awesomesauce_df

    # cooked_df = pd.read_csv("alignments.tsv", sep="\t")
    # cooked_df = cooked_df[
    #     (
    #         (cooked_df["query_chrom"] != cooked_df["break_chrom1"])
    #         & (cooked_df["query_chrom"] != cooked_df["break_chrom2"])
    #     )
    #     | (
    #         (
    #             (cooked_df["break_pos1"] > cooked_df["query_end"])
    #             | (cooked_df["break_pos1"] < cooked_df["query_pos"])
    #         )
    #         & (
    #             (cooked_df["break_pos2"] > cooked_df["query_end"])
    #             | (cooked_df["break_pos2"] < cooked_df["query_pos"])
    #         )
    #     )
    # ]
    # cooked_df = cooked_df.loc[
    #     cooked_df["query_cigar"].str.contains("S")
    #     | cooked_df["query_cigar"].str.contains("H")
    # ].drop(
    #     ["break_sv_type", "break_read_support", "break_features", "query_aln_full"],
    #     axis=1,
    # )
    # cooked_df.to_csv("minor_errors.tsv", sep="\t")

    # def process_row(row):
    #     if (
    #         row["break_pos1"] == row["query_pos"]
    #         or row["break_pos1"] == row["query_end"]
    #     ):
    #         if not row["break_pos1"] == row["query_pos"]:
    #             complement = {"C": "G", "G": "C", "A": "T", "T": "A"}
    #             return ("".join([complement[bp] for bp in row["query_aln_sub"]]))[::-1]
    #         else:
    #             return row["query_aln_sub"]
    #     else:
    #         if not row["break_pos2"] == row["query_end"]:
    #             complement = {"C": "G", "G": "C", "A": "T", "T": "A"}
    #             return ("".join([complement[bp] for bp in row["query_aln_sub"]]))[::-1]
    #         else:
    #             return row["query_aln_sub"]

    # new_awesomesauce_df["new_awesomesauce_query"] = new_awesomesauce_df.apply(
    #     process_row, axis=1
    # )
    # new_awesomesauce_df

    # def mini_homology_inator(str_list, dir):
    #     str_list = [i[::dir] for i in str_list]
    #     slice = 0
    #     for i in range(len(max(str_list, key=len))):
    #         temp_list = [i[:slice] for i in str_list]
    #         if not len(set(temp_list)) <= 1:
    #             break
    #         slice += 1
    #     return str_list[0][: slice - 1][::dir]

    # def homology_inator(first, last):
    #     first_str = mini_homology_inator(first, 1)
    #     last_str = mini_homology_inator(last, -1)
    #     longest_str = ""
    #     for i in range(min(len(first_str), len(last_str))):
    #         if first_str[:i] == last_str[-i:]:
    #             longest_str = first_str[:i]
    #     print(longest_str)
    #     return longest_str

    # homology_inator(first, last)

    # new_awesomesauce_df["first"] = new_awesomesauce_df.apply(lambda row: row["break_pos1"] == row["query_pos"] or row["break_pos1"] == row["query_end"], axis=1)
    # grouped_queries = new_awesomesauce_df.groupby(["break_pos1", "first"])["new_awesomesauce_query"].agg(list)
    # grouped_queries = grouped_queries.reset_index().pivot(index="break_pos1", columns="first", values="new_awesomesauce_query")
    # grouped_queries.columns = ["last", "first"]
    # grouped_queries['homology'] = grouped_queries.apply(lambda row: homology_inator(row["first"], row["last"]), axis=1)
    # grouped_queries
    # print(grouped_queries)
    #
    # grouped_queries['combined'] = grouped_queries['first'] + grouped_queries['last']
    # temp = grouped_queries[['combined', 'homology']].explode('combined')
    # homologies = dict(zip(temp['combined'], temp['homology']))
    # new_awesomesauce_df['homology'] = new_awesomesauce_df['new_awesomesauce_query'].apply(lambda x: homologies[x])
    # new_awesomesauce_df = new_awesomesauce_df[['break_chrom1', 'break_pos1', 'break_chrom2', 'break_pos2', 'AA_homology_seq', 'homology', 'query_name', 'split','query_chrom', 'query_pos', 'query_end', 'query_orientation', 'query_cigar', 'query_aln_sub', 'query_aln_full']]
    # new_awesomesauce_df = new_awesomesauce_df.reset_index(drop=True)
    # new_awesomesauce_df.to_csv('final.tsv', sep='\t')
    # new_awesomesauce_df

    # new_awesomesauce_df["homology_cut"] = new_awesomesauce_df.apply(lambda x: x['query_aln_sub'][:7] if x["break_start"] else x['query_aln_sub'][-7:], axis=1)
    # new_awesomesauce_df["wrong_dog"] = new_awesomesauce_df.apply(lambda x: str(x['AA_homology_seq']) not in x['homology_cut'], axis=1)
    # new_awesomesauce_df.loc[new_awesomesauce_df["wrong_dog"] == True]
