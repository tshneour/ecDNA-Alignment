import pandas as pd
import re
import argparse
import warnings
import pysam
from Bio.Seq import Seq
from Bio import SeqIO

warnings.filterwarnings("ignore")
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Refine AA SVs.")
    parser.add_argument(
        "file", type=str, help="alignments file to refine"
    )
    parser.add_argument(
        "fasta", type=str, help="fasta file used to generate bamfile"
    )
    parser.add_argument(
        "-l", "--list", help="shows all breakpoints", action="store_true"
    )
    parser.add_argument(
        "-b",
        "--breakpoints",
        type=int,
        help="list info for these breakpoints (use --list to see breakpoint_number)",
        nargs="+",
    )
    parser.add_argument(
        "-v", "--verbose", help="show information about specific reads", type=int
    )
    args = parser.parse_args()

    all_reads = pd.read_csv(args.file, sep="\t")
    fasta = pysam.FastaFile(args.fasta)

    # ------------- Define Functions ------------------

    def refine_step1(reads, certainty=1):
        beg_reg = r"^\d+(S|H)"
        end_reg = r"\d+(S|H)$"
        filtered_reads = reads[
            (reads["query_cigar"].str.contains("S"))
            | (reads["query_cigar"].str.contains("H"))
        ]

        filtered_reads["begin"] = filtered_reads.query_cigar.apply(
            lambda x: bool(re.search(beg_reg, x))
        )
        filtered_reads["end"] = filtered_reads.query_cigar.apply(
            lambda x: bool(re.search(end_reg, x))
        )

        filtered_reads["sv_end"] = filtered_reads.apply(
            lambda row: row["query_pos"]
            if row["begin"]
            else row["query_end"]
            if row["end"]
            else pd.NA,
            axis=1,
        )

        filtered_reads = filtered_reads.dropna(subset=["sv_end"])

        if args.verbose == 2:
            print(
                filtered_reads.sort_values(by="sv_end")[
                    [
                        "sv_end",
                        "query_name",
                        "proper_pair",
                        "query_pos",
                        "query_end",
                        "query_cigar",
                        "split",
                        "query_aln_full",
                    ]
                ].to_string(),
                "\n",
            )

        end_cands = filtered_reads["sv_end"].value_counts()
        best = end_cands.index[0] if end_cands.shape[0] > 0 else None

        if args.verbose == 1:
            verb_1 = filtered_reads[filtered_reads["sv_end"] == best][
                [
                    "sv_end",
                    "query_name",
                    "proper_pair",
                    "query_pos",
                    "query_end",
                    "query_cigar",
                    "split",
                    "query_aln_full",
                ]
            ]
            disc = verb_1[verb_1["proper_pair"] == "Discordant"]
            mates = all_reads[
                all_reads["query_name"].isin(disc.query_name)
                & ~all_reads["query_aln_full"].isin(disc.query_aln_full)
            ][
                [
                    "query_name",
                    "proper_pair",
                    "query_pos",
                    "query_end",
                    "query_cigar",
                    "split",
                    "query_aln_full",
                ]
            ]
            print(verb_1.to_string(), "\n")
            print()
            print("Mates of Discordant Reads from above Table: ")
            print(mates.to_string(), "\n")
        
        filtered_reads.split = filtered_reads.split.astype(int)
        group_rankings = filtered_reads.groupby("sv_end")['split'].agg(sum).sort_values(ascending=False)

        return end_cands, best, filtered_reads, group_rankings

    def refine_step2(best):
        def contains(read, break_cand):
            return break_cand in range(read["query_pos"] + 1, read["query_end"])

        # for break_cand in cands.index:
        #     display(break_cand, all_reads[all_reads.apply(lambda x: contains(x, break_cand), axis=1)])
        thing = all_reads.loc[all_reads.apply(lambda x: contains(x, best), axis=1)]
        if args.verbose:
            print(
                thing[
                    [
                        "query_name",
                        "query_pos",
                        "query_end",
                        "query_cigar",
                        "split",
                        "query_aln_full",
                    ]
                ].to_string()
                if len(thing) != 0
                else "None",
                "\n",
            )
        return len(thing)
    
    def check_overlap(left, left_rankings, right, right_rankings):
        left = left.reset_index(drop=True)
        right = right.reset_index(drop=True)

        left["split"] = left["split"].astype(bool)
        right["split"] = right["split"].astype(bool)
        left["query_short"] = left["query_short"].astype(str)
        right["query_short"] = right["query_short"].astype(str)

        # Eliminate unhelpful split reads
        left = left[(~left["split"].astype(bool)) | 
            (left["split"].astype(bool) & left["query_name"].isin(right.loc[right["split"].astype(bool), "query_name"]))]
        right = right[(~right["split"].astype(bool)) | 
            (right["split"].astype(bool) & right["query_name"].isin(left.loc[left["split"].astype(bool), "query_name"]))]

        # print(type(left))
        # print(left.columns)
        left["clip_len"] = left.apply(lambda x: len(x["query_aln_full"]) - len(x["query_aln_sub"]), axis=1)
        left["clipped"] = left.apply(lambda x: x["query_aln_full"][0:x["clip_len"]] if x["begin"] else x["query_aln_full"][-x["clip_len"] or len(x["query_aln_full"]):], axis=1)
        right["clip_len"] = right.apply(lambda x: len(x["query_aln_full"]) - len(x["query_aln_sub"]), axis=1)
        right["clipped"] = right.apply(lambda x: x["query_aln_full"][0:x["clip_len"]] if x["begin"] else x["query_aln_full"][-x["clip_len"] or len(x["query_aln_full"]):], axis=1)
        left_groups = left.groupby("sv_end")
        right_groups = right.groupby("sv_end")

        best_left = None
        best_right = None
        best_homology = None
        best_left_ins = None
        best_right_ins = None
        best_insertion = None

        for left_group, left_df in left_groups:
            for right_group, right_df in right_groups:
                first_row = left_df.loc[left_df["query_aln_sub"].str.len().idxmax()]
                last_row = right_df.loc[right_df["query_aln_sub"].str.len().idxmax()]
                first = first_row["query_aln_sub"] if first_row["end"] else reverse_complement(first_row["query_aln_sub"])
                last = last_row["query_aln_sub"] if last_row["begin"] else reverse_complement(last_row["query_aln_sub"])

                # Check for homology
                homology = homology_inator(first, last)
                # print(homology)
                hom_len = len(homology)
                first_nohomo = first[:-hom_len]
                last_nohomo = last[hom_len:]
                left_df["new_clipped"] = left_df.apply(lambda x: x["clipped"] if x["end"] else reverse_complement(x["clipped"]), axis=1)
                left_df["does_clip_match"] = left_df.apply(lambda x: len(homology_inator(x["new_clipped"], last_nohomo)) == len(x["new_clipped"]), axis=1)
                left_df.loc[left_df["split"], "does_clip_match"] = left_df.loc[left_df["split"]].apply(lambda x: x["query_short"] in (right_df.loc[right_df["split"], "query_short"]).values if 'H' in x["query_cigar"] else x["does_clip_match"], axis=1)
                right_df["new_clipped"] = right_df.apply(lambda x: x["clipped"] if x["begin"] else reverse_complement(x["clipped"]), axis=1)
                right_df["does_clip_match"] = right_df.apply(lambda x: len(homology_inator(first_nohomo, x["new_clipped"])) == len(x["new_clipped"]), axis=1)
                right_df.loc[right_df["split"], "does_clip_match"] = right_df.loc[right_df["split"]].apply(lambda x: x["query_short"] in (left_df.loc[left_df["split"], "query_short"]).values if 'H' in x["query_cigar"] else x["does_clip_match"], axis=1)
                # print(right_df["does_clip_match"].to_string())
                # print(left_df["new_clipped"].to_string())
                # print(last_nohomo)

                hom_sum_left = left_df["does_clip_match"].sum()
                hom_sum_right = right_df["does_clip_match"].sum()
                # print(sum_left, sum_right)

                # Check for insertion
                first_row = left_df.loc[left_df["new_clipped"].str.len().idxmax()]
                last_row = right_df.loc[right_df["new_clipped"].str.len().idxmax()]
                first_clipped = first_row["new_clipped"]
                last_clipped = last_row["new_clipped"]
                insertion = homology_inator(last_clipped, first_clipped)
                insertion_len = len(insertion)
                left_df["does_clip_match"] = left_df.apply(lambda x: len(homology_inator(x["new_clipped"][insertion_len:], last)) == len(x["new_clipped"][insertion_len:]), axis=1)
                left_df.loc[left_df["split"], "does_clip_match"] = left_df.loc[left_df["split"]].apply(lambda x: x["query_short"] in (right_df.loc[right_df["split"], "query_short"]).values if 'H' in x["query_cigar"] else x["does_clip_match"], axis=1)
                right_df["does_clip_match"] = right_df.apply(lambda x: len(homology_inator(first, x["new_clipped"][:-insertion_len])) == len(x["new_clipped"][:-insertion_len]), axis=1)
                right_df.loc[right_df["split"], "does_clip_match"] = right_df.loc[right_df["split"]].apply(lambda x: x["query_short"] in (left_df.loc[left_df["split"], "query_short"]).values if 'H' in x["query_cigar"] else x["does_clip_match"], axis=1)
                ins_sum_left = left_df["does_clip_match"].sum()
                ins_sum_right = right_df["does_clip_match"].sum()

                if (hom_sum_left == len(left_df) and hom_sum_right == len(right_df)):
                    print("Yes! All soft clips match!")
                    if (best_left is not None): 
                        print("Hmmmm.... Multiple candidates with matching soft clips", f"Namely left: {left_group} right: {right_group} homology: {homology}")
                    best_left = left_group
                    best_right = right_group
                    best_homology = homology if first_row["end"] else reverse_complement(homology)
                elif (hom_sum_left > 0 or hom_sum_right > 0):
                    print("Hmmmmm...... only some soft clips matched", f"Namely left: {hom_sum_left}/{len(left_df)} right: {hom_sum_right}/{len(right_df)}")
                    print("More info:", f"Namely left: {left_group} right: {right_group} homology: {homology}")
                    print(left_df[["query_short", "query_chrom", "query_pos", "query_cigar", "query_aln_full", "query_aln_sub"]].to_string())
                    print(right_df[["query_short", "query_chrom", "query_pos", "query_cigar", "query_aln_full", "query_aln_sub"]].to_string())
                    print()


                if (insertion_len == 0):
                    continue
                # print(first, last, "homology:", homology)
                # print(right_df[["does_clip_match", "query_cigar", "new_clipped"]].to_string())
                if (ins_sum_left == len(left_df) and ins_sum_right == len(right_df)):
                    print()
                    print("Trying as insertion")
                    print("Yes! All soft clips match!")
                    # print("INSERTION:", insertion)
                    if (best_left is not None): 
                        print("Hmmmm.... Multiple candidates with matching soft clips", f"Namely left: {left_group} right: {right_group} insertion: {insertion}")
                    best_left_ins = left_group
                    best_right_ins = right_group
                    best_insertion = insertion if first_row["end"] else reverse_complement(insertion)
                elif (ins_sum_left > 0 or ins_sum_right > 0):
                    print()
                    print("Trying as insertion")
                    print("Hmmmmm...... only some soft clips matched", f"Namely left: {ins_sum_left}/{len(left_df)} right: {ins_sum_right}/{len(right_df)}")
                    print("More info:", f"Namely left: {left_group} right: {right_group} insertion: {insertion}")
                    print(left_df[["query_short", "query_chrom", "query_pos", "query_cigar", "query_aln_full", "query_aln_sub"]].to_string())
                    print(right_df[["query_short", "query_chrom", "query_pos", "query_cigar", "query_aln_full", "query_aln_sub"]].to_string())
                    print()
        
        return best_left, best_right, best_homology, best_left_ins, best_right_ins, best_insertion
        
    def reverse_complement(dna: str) -> str:
        """Returns the reverse complement of a DNA sequence using Biopython."""
        return str(Seq(dna).reverse_complement())

    def homology_inator(first, last):
        first_str = first
        last_str = last
        longest_str = ""
        for i in range(min(len(first_str), len(last_str)) + 1):
            if first_str[-i:] == last_str[:i]:
                longest_str = first_str[-i:]
        # print(first, last, longest_str)
        return longest_str
    
    def extract_sequence(seq_id, start, end):

        # Check if the sequence ID exists in the FAI index
        if seq_id in fasta.references:
            # Fetch the sequence from the FASTA file in the specified range (start, end)
            sequence = fasta.fetch(seq_id, start, end)
            return sequence
        else:
            print("UHH OHH AAAAAAAAAAAA")
            return None  # If the sequence ID is not found
    
   # ------------- Do refinement ------------------

    # all_reads.loc[all_reads['query_cigar'].str.contains('H'), 'query_aln_full'] = all_reads.loc[all_reads['query_cigar'].str.contains('H')].apply(
    #     lambda row: extract_sequence(row["query_chrom"], row['query_pos'], row['query_end']), axis=1
    # )

    pd.set_option("mode.chained_assignment", None)
    svs = all_reads.groupby("break_pos1")
    for group_name, sv in svs:
        left = sv[
            (sv["query_pos"] >= (sv["break_pos1"] - 350))
            & (sv["query_end"] <= (sv["break_pos1"] + 350))
        ]
        right = sv[
            (sv["query_pos"] >= (sv["break_pos2"] - 350))
            & (sv["query_end"] <= (sv["break_pos2"] + 350))
        ]
        if args.list:
            out = (
                svs[
                    [
                        "break_chrom1",
                        "break_pos1",
                        "break_chrom2",
                        "break_pos2",
                        "AA_homology_seq",
                    ]
                ]
                .first()
                .reset_index(drop=True)
            )
            out.index.name = "breakpoint_index"
            print(out.to_string())
            break
        if args.breakpoints:
            breaks = (
                svs["break_pos1"]
                .first()
                .reset_index(drop=True)
                .iloc[args.breakpoints]
                .to_list()
            )
            if group_name not in breaks:
                continue

        if args.verbose:
            print(
                f"Reads for {'best' if args.verbose == 1 else ('' if args.verbose == 2 else '')} left SV candidates"
            )
        left_cands, best_left, left_groups, left_rankings = refine_step1(left)
        if args.verbose:
            print(
                f"Reads for {'best' if args.verbose == 1 else ('' if args.verbose == 2 else '')} right SV candidates"
            )
        right_cands, best_right, right_groups, right_rankings = refine_step1(right)

        new_best_left, new_best_right, homology, new_best_insL, new_best_insR, insertion = check_overlap(left_groups, left_rankings, right_groups, right_rankings)
        print(f"New refinement (Homology) -> left: {new_best_left} right: {new_best_right} homology: {homology}")
        print(f"New refinement (Insertion) -> left: {new_best_insL} right: {new_best_insR} insertion: {insertion}")

        # print(left_groups.to_string())
        # print(right_groups.to_string())

        print("Original AA breakpoint:")
        print(
            sv.iloc[0][
                [
                    "break_chrom1",
                    "break_pos1",
                    "break_chrom2",
                    "break_pos2",
                    "AA_homology_seq",
                    "break_orientation"
                ]
            ].to_string(),
            "\n",
        )

        print("Best left candidate:", best_left)
        print("Best right candidate:", best_right, "\n")

        if args.verbose:
            print("Overlapping reads with best left candidate:")
        left_ovrs = refine_step2(best_left)
        if args.verbose:
            print("Overlapping reads with best right candidate:")
        right_ovrs = refine_step2(best_right)

        print("Reads overlapping on best left candidate:", left_ovrs)
        print("Reads overlapping on best right candidate:", right_ovrs)
        print()