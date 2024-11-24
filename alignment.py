import pandas as pd
import pysam
import re
import argparse

def get_pairs(regions, chrom1, chrom2):
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
            if any(aln in tup for tup in split_alignments) or any(aln in tup for tup in nonsplit_alignments):
                print("Duplicate found")
                continue

            mate = samfile.mate(aln)

            if aln.is_supplementary or aln.has_tag("SA"):
                sa_tag = (aln.get_tag("SA")).split(";")[:-1]
                sa_tag = sa_tag[0]
                sup_name, sup_pos, sup_strand, sup_cigar, sup_mapq, sup_tlen = sa_tag.split(",")
                sup_pos = int(sup_pos)
                sup_tlen = int(sup_tlen)

                found_split2 = False
                for sup_read in samfile.fetch(sup_name, sup_pos, sup_pos + 1):
                    if (sup_read.query_name == aln.query_name and sup_read.reference_name == sup_name and (sup_read.reference_start + 1) == sup_pos and re.sub(r'[SH]', '', sup_read.cigarstring) == re.sub(r'[SH]', '', sup_cigar)):
                        non_sup = aln if not aln.is_supplementary else sup_read
                        sup = aln if aln.is_supplementary else sup_read
                        read1 = non_sup if non_sup.is_read1 else mate
                        read2 = non_sup if non_sup.is_read2 else mate
                        split_alignments.append((read1, sup, read2))
                        found_split2 = True
                        break
                if not found_split2:
                    print(sup_name, sup_pos, sup_cigar, sup_tlen)
                    raise NameError("split mate not found")

            elif mate.is_supplementary or mate.has_tag("SA"):
                sa_tag = (mate.get_tag("SA")).split(";")[:-1]
                sa_tag = sa_tag[0]
                sup_name, sup_pos, sup_strand, sup_cigar, sup_mapq, sup_tlen = sa_tag.split(",")
                sup_pos = int(sup_pos)
                sup_tlen = int(sup_tlen)

                found_split2 = False
                for sup_read in samfile.fetch(sup_name, sup_pos, sup_pos + 1):
                    if (sup_read.query_name == mate.query_name and sup_read.reference_name == sup_name and (sup_read.reference_start + 1) == sup_pos and re.sub(r'[SH]', '', sup_read.cigarstring) == re.sub(r'[SH]', '', sup_cigar)):
                        non_sup = mate if not mate.is_supplementary else sup_read
                        sup = mate if mate.is_supplementary else sup_read
                        read1 = non_sup if non_sup.is_read1 else aln
                        read2 = non_sup if non_sup.is_read2 else aln
                        split_alignments.append((read1, sup, read2))
                        found_split2 = True
                        break
                if not found_split2:
                    print(sup_name, sup_pos, sup_cigar, sup_tlen)
                    raise NameError("split mate not found")
                
            elif not aln.is_proper_pair and aln.next_reference_name == (chrom2 if i == 0 else chrom1):
                read1 = aln if aln.is_read1 else mate
                read2 = aln if aln.is_read2 else mate
                nonsplit_alignments.append((read1, read2))
            
            elif aln.is_proper_pair:
                read1 = aln if aln.is_read1 else mate
                read2 = aln if aln.is_read2 else mate
                nonsplit_alignments.append((read1, read2))
    
    return split_alignments, nonsplit_alignments

def fetch_alignments(bamfile, chrom1, pos1, chrom2, pos2, sv_type, read_support, features, orientation, hom_len, hom, region_size=150, mapq_threshold=15):
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

    split_alignments, nonsplit_alignments = get_pairs([region1, region2], chrom1, chrom2)

    # Collect details of split alignments
    split_alignment_details = []
    for pair in split_alignments:
        read1_1, read1_2, read2 = pair

        # Filter out low quality reads
        if  ((read1_1.mapping_quality > mapq_threshold and read1_2.mapping_quality > mapq_threshold) or read2.mapping_quality > mapq_threshold) or ((read1_1.is_mapped and read1_2.is_mapped) or read2.is_mapped):
            for read in pair:
                split_alignment_details.append({
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
                    "query_short": read.query_name.split(':')[-1],
                    "split": read.has_tag("SA"),
                    "proper_pair": "Concordant" if read.is_proper_pair else "Discordant",
                    "read_num": "1" if read.is_read1 else "2",
                    "query_chrom": read.reference_name,
                    "query_pos": read.reference_start + 1,
                    "query_end": read.reference_end,
                    "query_orientation": "+" if read.is_forward else "-",
                    "query_cigar": read.cigarstring,
                    "query_aln_full": read.query_sequence,
                    "query_aln_sub": read.query_alignment_sequence
                })

    # Collect details of nonsplit alignments
    nonsplit_alignment_details = []
    for pair in nonsplit_alignments:
        read1, read2 = pair
        # Filter out low quality reads
        if  (read1.mapping_quality > mapq_threshold or read2.mapping_quality > mapq_threshold) or (read1.is_mapped or read2.is_mapped):
            for read in pair:
                nonsplit_alignment_details.append({
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
                    "query_short": read.query_name.split(':')[-1],
                    "split": read.has_tag("SA"),
                    "proper_pair": "Concordant" if read1.is_proper_pair else "Discordant",
                    "read_num": "1" if read.is_read1 else "2",
                    "query_chrom": read.reference_name,
                    "query_pos": read.reference_start + 1,
                    "query_end": read.reference_end,
                    "query_orientation": "+" if read.is_forward else "-",
                    "query_cigar": read.cigarstring,
                    "query_aln_full": read.query_sequence,
                    "query_aln_sub": read.query_alignment_sequence
                })
        
    samfile.close()

    # Convert details to pandas DataFrames
    split_df = pd.DataFrame(split_alignment_details)
    nonsplit_df = pd.DataFrame(nonsplit_alignment_details)
    # paired_df = pd.DataFrame(paired_alignment_details)

    return split_df, nonsplit_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collect reads from the ends of a SV.")
    parser.add_argument("--sum", help="path to AA amplicon summaries file", type=str, default="./K562/K562_summaries.tsv")
    parser.add_argument("--bam", help="path to bamfile", type=str, default="./K562/K562_hg19_cnvkit.cs.rmdup.bam")
    args = parser.parse_args()
    # Define the BAM file and positions
    bamfile = args.bam
    samfile = pysam.AlignmentFile(bamfile, "rb")
    df = pd.read_csv(args.sum, sep="\t")
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
        split_df, nonsplit_df = fetch_alignments(bamfile, chrom1, pos1, chrom2, pos2, sv_type, read_support, features, orientation, hom_len, homology, 150)
        num_split += len(split_df)/3
        num_paired += (2 * len(split_df) / 3) + len(nonsplit_df)
        # num_paired += len(paired_df)
        

        print(f"Reads at junction defined by {chrom1} {pos1} and {chrom2} {pos2}:")
        print(f"Number split pairs: {len(split_df)/3}")
        print(f"Number nonsplit pairs: {len(nonsplit_df)/2}")
        print(f"Number paired: {(2 * len(split_df) / 3) + len(nonsplit_df)}")
        output = pd.concat([output, split_df], ignore_index=True)
        output = pd.concat([output, nonsplit_df], ignore_index=True)

    output.to_csv('alignments.tsv', sep='\t', index=False)

    print("Total Number of split reads:", num_split)
    print("Total Number of paired reads:", num_paired)