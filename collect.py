#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam
import re
import argparse
import os

# import gzip
import subprocess

_sh_re = re.compile(r"[SH]")


# Identify sv orientation from split alignment cigar
def get_sv_ori(cigar):
    l_matches = re.search(r"^\d+[SH]", cigar)
    r_matches = re.search(r"\d+[SH]$", cigar)
    # print(cigar)
    # print(l_matches)
    # print(r_matches)

    if l_matches:
        return "-"
    elif r_matches:
        return "+"
    else:
        raise Exception("No clip matches in cigar string")


# Use to filter out non-supportive discordant and split pairs
def filter_for_pos_ori(
    read1,
    read2,
    chrom1,
    break_pos1,
    chrom2,
    break_pos2,
    sv_type,
    sv_orientation,
    sup=None,
):
    left, right = None, None

    # If split read, treat read1 as primary alignment and read2 as supplementary alignment
    if sup is not None:
        read1 = read1 if sup.is_read1 else read2
        read2 = sup

    # Check that split alignments/discordant mates are on right chromosomes
    if not (
        (read1.reference_name == chrom1 and read2.reference_name == chrom2)
        or (read1.reference_name == chrom2 and read2.reference_name == chrom1)
    ):
        return False

    # Determine which alignment represents left side of sv and which represents right side
    if chrom1 == chrom2:
        left = (
            read1
            if abs(read1.reference_start - break_pos1)
            < abs(read2.reference_start - break_pos1)
            else read2
        )
        right = (
            read1
            if abs(read1.reference_start - break_pos2)
            < abs(read2.reference_start - break_pos2)
            else read2
        )
    else:
        left = read1 if read1.reference_name == chrom1 else read2
        right = read1 if read1.reference_name == chrom2 else read2

    # Check that splits/mates are reasonably close to AA SV ends
    if (
        abs(left.reference_start - break_pos1) >= 500
        or abs(right.reference_start - break_pos2) >= 500
    ):
        return False

    # Left inversion
    if sv_orientation == "--":
        # Expected AA SV direction/Illumina orientations
        if sup is not None and not (
            get_sv_ori(left.cigarstring) == "-" and get_sv_ori(right.cigarstring) == "-"
        ):
            return False
        if (
            sup is None
            and sv_type != "interchromosomal"
            and not (left.is_reverse and right.is_reverse)
        ):
            return False

    # Right inversion
    elif sv_orientation == "++":
        # Expected AA SV direction/Illumina orientations
        if sup is not None and not (
            get_sv_ori(left.cigarstring) == "+" and get_sv_ori(right.cigarstring) == "+"
        ):
            return False
        if (
            sup is None
            and sv_type != "interchromosomal"
            and not (left.is_forward and right.is_forward)
        ):
            return False
    # Deletion
    elif sv_orientation == "+-":
        # Expected AA SV direction/Illumina orientations
        if sup is not None and not (
            get_sv_ori(left.cigarstring) == "+" and get_sv_ori(right.cigarstring) == "-"
        ):
            return False
        if (
            sup is None
            and sv_type != "interchromosomal"
            and not (left.is_forward and right.is_reverse)
        ):
            return False
    # Duplication
    elif sv_orientation == "-+":
        # Expected AA SV direction/Illumina orientations
        if sup is not None and not (
            get_sv_ori(left.cigarstring) == "-" and get_sv_ori(right.cigarstring) == "+"
        ):
            return False
        if (
            sup is None
            and sv_type != "interchromosomal"
            and not (left.is_reverse and right.is_forward)
        ):
            return False

    return True


def check_strict(read1, read2, pos1, pos2, region_size, sup=None):
    def overlaps_region(read, center, region_size):
        region_start = center - region_size
        region_end = center + region_size
        return read.reference_start <= region_end and read.reference_end >= region_start

    def overlaps_either_region(read):
        return overlaps_region(read, pos1, region_size) or overlaps_region(
            read, pos2, region_size
        )

    # Check reads overlap with one of either region
    if sup is not None:
        return (
            overlaps_either_region(read1)
            and overlaps_either_region(read2)
            and overlaps_either_region(sup)
        )
    else:
        return overlaps_either_region(read1) and overlaps_either_region(read2)


def rev_comp(dna):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    reversed_dna = dna[::-1]
    return "".join(complement[base] for base in reversed_dna)


def get_pairs(
    regions, chrom1, chrom2, break_pos1, break_pos2, sv_type, sv_orientation, args
):
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
    split_leftover_alignments = []

    # Collect alignments from left and right regions
    reads_by_name = {}  # query_name -> list[AlignedSegment]
    for region in regions:
        # Remove optical duplicates and unmapped reads
        for aln in region:
            if aln.is_duplicate or not aln.is_mapped:
                continue
            reads_by_name.setdefault(aln.query_name, []).append(aln)

    def passes_pos_filter(read1, read2, sup=None):
        return filter_for_pos_ori(
            read1,
            read2,
            chrom1,
            break_pos1,
            chrom2,
            break_pos2,
            sv_type,
            sv_orientation,
            sup,
        )

    # 3) iterate each group (all alignments for the same query_name) and resolve pairs locally
    for qname, alns in reads_by_name.items():
        # Remove duplicates
        alns_set = set(alns)
        # find primary alignments and any supplementary/secondary ones
        primary = [a for a in alns_set if not (a.is_supplementary or a.is_secondary)]
        sup_or_sec = [a for a in alns_set if (a.is_supplementary or a.is_secondary)]

        # If we have complete pairs, assign to lists
        if len(primary) == 2:
            # Assign the nonsplit pairs
            a = primary[0]
            b = primary[1]
            read1 = a if a.is_read1 else b if b.is_read1 else a
            read2 = b if b.is_read2 else a if a.is_read2 else b

            if (
                not read1.has_tag("SA")
                and not read2.has_tag("SA")
                and not sup_or_sec
                and read1.is_proper_pair
            ):
                nonsplit_alignments.append((read1, read2))
                continue

            if (
                not read1.has_tag("SA")
                and not read2.has_tag("SA")
                and not sup_or_sec
                and not read1.is_proper_pair
                and passes_pos_filter(read1, read2)
            ):
                nonsplit_alignments.append((read1, read2))
                nonsplit_alignments.append((read1, read2))
                continue
            # Assign pairs with split
            if primary and sup_or_sec:
                if len(sup_or_sec) == 1:
                    prim = read1 if read1.has_tag("SA") else read2
                    sup = sup_or_sec[0]

                    if sup.query_name != prim.query_name:
                        print(
                            "Mismatch between sup and prim alignments. This shouldn't happen..."
                        )
                        continue

                    if passes_pos_filter(read1, read2, sup):
                        split_alignments.append((read1, sup, read2))
                        split_alignments.append((read1, sup, read2))
                        continue

                elif len(sup_or_sec) == 2:
                    sup1, sup2 = (
                        sup_or_sec[0],
                        sup_or_sec[1]
                    ) if sup_or_sec[0].is_read1 else (sup_or_sec[1], sup_or_sec[0])

                    if sup1.is_read1 == sup2.is_read1:
                        split_leftover_alignments.append((read1, sup1, sup2, read2) if sup1.is_read1 else (read1, read2, sup1, sup2))
                        continue

                    if (
                        sup1.query_name != read1.query_name
                        or sup2.query_name != read2.query_name
                    ):
                        print(
                            "Mismatch between sup and prim alignments. This shouldn't happen..."
                        )
                        print(read1.query_name, sup1.query_name, read2.query_name, sup2.query_name)
                        continue

                    if passes_pos_filter(read1, read2, sup1) and passes_pos_filter(
                        read1, read2, sup2
                    ):
                        split_alignments.append((read1, sup1, read2, sup2))
                        split_alignments.append((read1, sup1, read2, sup2))
                        continue
                else:
                    sups1 = [sup for sup in sup_or_sec if sup.is_read1]
                    sups2 = [sup for sup in sup_or_sec if sup.is_read2]
                    result = tuple([read1] + sups1 + [read2] + sups2)
                    split_leftover_alignments.append(result)

    return split_alignments, nonsplit_alignments, split_leftover_alignments


# @profile
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
    args,
    sample,
    samfile,
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
    - args.refine: int, the size of the region around the positions to fetch alignments from.
    - mapq_threshold: int, cutoff for collecting reads above certain mapping quality.

    Returns:
    - split_df: DataFrame, details of split alignments.
    - paired_df: DataFrame, details of discordant alignments.
    """

    # Fetch alignments args.refine away from both positions
    region1 = list(samfile.fetch(chrom1, pos1 - args.refine, pos1 + args.refine + 1))
    region2 = list(samfile.fetch(chrom2, pos2 - args.refine, pos2 + args.refine + 1))

    split_alignments, nonsplit_alignments, split_leftover_alignments = get_pairs(
        [region1, region2], chrom1, chrom2, pos1, pos2, sv_type, orientation, args
    )

    # Collect details of split alignments
    split_alignment_details = []
    for pair in split_alignments:
        if len(pair) == 3:
            read1_1, read1_2, read2 = pair
            if args.strict and not check_strict(
                read1_1, read2, pos1, pos2, args.refine, sup=read1_2
            ):
                continue
            # Filter out low quality reads
            if (
                (
                    read1_1.mapping_quality > mapq_threshold
                    and read1_2.mapping_quality > mapq_threshold
                )
                or read2.mapping_quality > mapq_threshold
            ) or ((read1_1.is_mapped and read1_2.is_mapped) or read2.is_mapped):
                for read in pair:
                    query_aln_full = read.query_sequence
                    if read.is_secondary or read.is_supplementary:
                        prim_read = read1_1 if read1_1.has_tag("SA") else read2
                        query_aln_full = (
                            prim_read.query_sequence
                            if prim_read.is_forward == read.is_forward
                            else rev_comp(prim_read.query_sequence)
                        )

                    record = {
                        "break_chrom1": chrom1,
                        "break_pos1": pos1 + 1,
                        "break_chrom2": chrom2,
                        "break_pos2": pos2 + 1,
                        "break_sv_type": sv_type,
                        "break_orientation": orientation,
                        "homology_len": hom_len,
                        "homology_seq": hom,
                        "query_name": read.query_name,
                        "query_short": ":".join(read.query_name.split(":")[-2:]),
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
                        "query_aln_full": query_aln_full,
                        "query_aln_sub": read.query_alignment_sequence,
                        "sample": sample,
                        "query_qualities": read.query_qualities,
                    }
                    if read_support is not None:
                        record["break_read_support"] = read_support
                    if features is not None:
                        record["break_features"] = features
                    split_alignment_details.append(record)
        if len(pair) == 4:
            read1_1, read1_2, read2_1, read2_2 = pair
            if args.strict and (
                not check_strict(read1_1, read2_1, pos1, pos2, args.refine, sup=read1_2)
                or not check_strict(
                    read1_1, read2_1, pos1, pos2, args.refine, sup=read2_2
                )
            ):
                continue
            # Filter out low quality reads
            if (
                read1_1.mapping_quality > mapq_threshold
                and read1_2.mapping_quality > mapq_threshold
                and read2_1.mapping_quality > mapq_threshold
                and read2_2.mapping_quality > mapq_threshold
            ) or (
                (read1_1.is_mapped and read1_2.is_mapped)
                or (read2_1.is_mapped and read2_2.is_mapped)
            ):
                for read in pair:
                    query_aln_full = read.query_sequence
                    if read.is_secondary or read.is_supplementary:
                        prim_read = (
                            read1_1 if (read1_1.is_read1 == read.is_read1) else read2_1
                        )
                        query_aln_full = (
                            prim_read.query_sequence
                            if prim_read.is_forward == read.is_forward
                            else rev_comp(prim_read.query_sequence)
                        )

                    record = {
                        "break_chrom1": chrom1,
                        "break_pos1": pos1 + 1,
                        "break_chrom2": chrom2,
                        "break_pos2": pos2 + 1,
                        "break_sv_type": sv_type,
                        "break_orientation": orientation,
                        "homology_len": hom_len,
                        "homology_seq": hom,
                        "query_name": read.query_name,
                        "query_short": ":".join(read.query_name.split(":")[-2:]),
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
                        "query_aln_full": query_aln_full,
                        "query_aln_sub": read.query_alignment_sequence,
                        "sample": sample,
                        "query_qualities": read.query_qualities,
                    }
                    if read_support is not None:
                        record["break_read_support"] = read_support
                    if features is not None:
                        record["break_features"] = features
                    split_alignment_details.append(record)
    
    split_leftover_alignment_details = []
    for pair in split_leftover_alignments:
        reads = list(pair)
        read1all = [read for read in reads if read.is_read1]
        read2all = [read for read in reads if read.is_read2]
        read1 = read1all[0]
        read2 = read2all[0]
        read1allsorted = sorted(read1all, key= lambda read: last_match_coord(read.cigartuples[::-1] if read1.is_forward != read.is_forward else read.cigartuples))
        read2allsorted = sorted(read2all, key= lambda read: last_match_coord(read.cigartuples[::-1] if read2.is_forward != read.is_forward else read.cigartuples))
        # Filter out low quality reads
        if all(read.is_mapped and read.mapping_quality > mapq_threshold for read in reads):
            for read in read1allsorted + read2allsorted:
                query_aln_full = read.query_sequence
                if read.is_secondary or read.is_supplementary:
                    prim_read = (
                        read1 if read.is_read1 else read2
                    )
                    query_aln_full = (
                        prim_read.query_sequence
                        if prim_read.is_forward == read.is_forward
                        else rev_comp(prim_read.query_sequence)
                    )

                record = {
                    "break_chrom1": chrom1,
                    "break_pos1": pos1 + 1,
                    "break_chrom2": chrom2,
                    "break_pos2": pos2 + 1,
                    "break_sv_type": sv_type,
                    "break_orientation": orientation,
                    "homology_len": hom_len,
                    "homology_seq": hom,
                    "query_name": read.query_name,
                    "query_short": ":".join(read.query_name.split(":")[-2:]),
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
                    "query_aln_full": query_aln_full,
                    "query_aln_sub": read.query_alignment_sequence,
                    "sample": sample,
                    "query_qualities": read.query_qualities,
                }
                if read_support is not None:
                    record["break_read_support"] = read_support
                if features is not None:
                    record["break_features"] = features
                split_leftover_alignment_details.append(record)

    # Collect details of nonsplit alignments
    nonsplit_alignment_details = []
    for pair in nonsplit_alignments:
        read1, read2 = pair
        if args.strict and not check_strict(read1, read2, pos1, pos2, args.refine):
            continue
        # Filter out low quality reads
        if (
            read1.mapping_quality > mapq_threshold
            or read2.mapping_quality > mapq_threshold
        ) or (read1.is_mapped or read2.is_mapped):
            for read in pair:
                record = {
                    "break_chrom1": chrom1,
                    "break_pos1": pos1 + 1,
                    "break_chrom2": chrom2,
                    "break_pos2": pos2 + 1,
                    "break_sv_type": sv_type,
                    "break_orientation": orientation,
                    "homology_len": hom_len,
                    "homology_seq": hom,
                    "query_name": read.query_name,
                    "query_short": ":".join(read.query_name.split(":")[-2:]),
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
                    "sample": sample,
                    "query_qualities": read.query_qualities,
                }
                if read_support is not None:
                    record["break_read_support"] = read_support
                if features is not None:
                    record["break_features"] = features
                nonsplit_alignment_details.append(record)

    # Convert details to pandas DataFrames
    split_df = pd.DataFrame(split_alignment_details)
    nonsplit_df = pd.DataFrame(nonsplit_alignment_details)
    split_leftover_df = pd.DataFrame(split_leftover_alignment_details)

    return split_df, nonsplit_df, split_leftover_df

def last_match_coord(cigartuple):
    count = 0
    maxmatch = 0
    for pair in cigartuple:
        op, num = pair
        if op != pysam.CIGAR_OPS.CMATCH:
            count += num
        else:
            maxmatch = count
    return maxmatch

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
    parser.add_argument(
        "--strict",
        help="Force collection of read pairs which fully align within the region of interest",
        action="store_true",
    )
    parser.add_argument("-f", "--file", help="File to store alignments", type=str)
    args = parser.parse_args()
    # Define the BAM file and positions
    bamfile = args.bam
    samfile = pysam.AlignmentFile(bamfile, "rb")
    # print(samfile.header)
    df = pd.concat(
        [
            pd.read_csv(os.path.join(args.sum, f), sep="\t").assign(
                sample=f.split("_")[1] if "amplicon" in f else f
            )
            for f in os.listdir(args.sum)
            if f.endswith(".tsv")
        ],
        ignore_index=True,
    )

    # Track which optional columns are present in the input
    has_read_support = "read_support" in df.columns
    has_features = "features" in df.columns

    required_cols = ["chrom1", "pos1", "chrom2", "pos2", "sv_type", "orientation"]
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        raise RuntimeError(f"Missing required columns: {missing}")

    has_homology = (
        "homology_length" in df.columns and
        "homology_sequence" in df.columns
    )

    # Define output table
    output = pd.DataFrame()
    num_split = 0
    num_discord = 0
    num_concord = 0
    num_paired = 0

    output = []
    leftover_output = []
    for index, row in enumerate(df.itertuples(index=False), 1):
        print(f"Breakpoint {index}")
        chrom1 = row.chrom1
        pos1 = row.pos1
        chrom2 = row.chrom2
        pos2 = row.pos2
        sv_type = row.sv_type
        read_support = row.read_support if has_read_support else None
        features = row.features if has_features else None
        orientation = row.orientation
        if has_homology:
            hom_len = row.homology_length
            homology = row.homology_sequence
        else:
            hom_len = None
            homology = None
        sample = row.sample

        # Fetch the DataFrames for split alignments and paired alignments
        split_df, nonsplit_df, split_leftover_df = fetch_alignments(
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
            args,
            sample,
            samfile,
        )

        if args.verbose:
            split_len = len(split_df)
            nonsplit_len = len(nonsplit_df)
            num_split += split_len / 3
            num_paired += (2 * split_len / 3) + nonsplit_len
            print(f"Reads at junction defined by {chrom1} {pos1} and {chrom2} {pos2}:")
            print(f"Number split pairs: {split_len / 3}")
            print(f"Number nonsplit pairs: {nonsplit_len / 2}")
            print(f"Number paired: {(2 * split_len / 3) + nonsplit_len}")
        if not split_df.empty:
            output.append(split_df)
        if not nonsplit_df.empty:
            output.append(nonsplit_df)
        if not split_leftover_df.empty:
            leftover_output.append(split_leftover_df)

    samfile.close()

    output = pd.concat(output, ignore_index=True) if output else pd.DataFrame()
    leftover_output = pd.concat(leftover_output, ignore_index=True) if leftover_output else pd.DataFrame()

    # Get split reads to determine how to filter concordant reads
    output_copy = output.copy()
    mask_split = output_copy["split"]
    output_copy = output_copy.loc[mask_split]

    mask_sh = output_copy["query_cigar"].str.endswith(("S", "H"))
    output_copy["ref_pos"] = np.where(
        mask_sh, output_copy["query_end"], output_copy["query_pos"]
    )

    # Determine sv end (left or right) of split alignments
    left_is_first = (
        (output_copy["break_chrom1"] != output_copy["break_chrom2"])
        & (output_copy["query_chrom"] == output_copy["break_chrom1"])
    ) | (
        (output_copy["break_chrom1"] == output_copy["break_chrom2"])
        & (
            (output_copy["query_pos"] - output_copy["break_pos1"]).abs()
            < (output_copy["query_pos"] - output_copy["break_pos2"]).abs()
        )
    )
    output_copy["is_left"] = left_is_first

    sv_cols = ["break_chrom1", "break_pos1"]

    # Turn groups of split reads into aggregated series
    def sv_summary(grp):
        mask_left = grp["is_left"]
        mask_right = ~mask_left

        def most_common(series):
            if series is None or series.empty:
                return pd.NA
            vc = series.value_counts()
            if vc.empty:
                return pd.NA
            return vc.idxmax()

        left_pos = most_common(grp.loc[mask_left, "ref_pos"])
        right_pos = most_common(grp.loc[mask_right, "ref_pos"])

        all_left_within_10 = (
            grp.loc[mask_left, "ref_pos"]
            .sub(grp.loc[mask_left, "break_pos1"])
            .abs()
            .le(10)
            .all()
        )
        all_right_within_10 = (
            grp.loc[mask_right, "ref_pos"]
            .sub(grp.loc[mask_right, "break_pos2"])
            .abs()
            .le(10)
            .all()
        )

        return pd.Series(
            {
                "left_pos": left_pos,
                "right_pos": right_pos,
                "all_left_within_10": all_left_within_10,
                "all_right_within_10": all_right_within_10,
            }
        )

    if output_copy.shape[0] != 0:
        df_sv = output_copy.groupby(sv_cols).apply(sv_summary).reset_index()

        df_sv = df_sv[(df_sv["all_left_within_10"]) & (df_sv["all_right_within_10"])].copy()
        out = output.merge(df_sv, on=sv_cols, how="left")

        # Determine sv end (left or right) of split alignments
        left_is_first = (
            (out["break_chrom1"] != out["break_chrom2"])
            & (out["query_chrom"] == out["break_chrom1"])
        ) | (
            (out["break_chrom1"] == out["break_chrom2"])
            & (
                (out["query_pos"] - out["break_pos1"]).abs()
                < (out["query_pos"] - out["break_pos2"]).abs()
            )
        )
        out["is_left"] = left_is_first

        drop_mask = pd.Series(False, index=out.index)

        is_concordant = out["proper_pair"] == "Concordant"
        ori0 = out["break_orientation"].str.get(0)
        ori1 = out["break_orientation"].str.get(1)

        left_present = out["left_pos"].notna()
        cond_left_plus = (
            is_concordant
            & left_present
            & out["is_left"]
            & (ori0 == "+")
            & (out["query_end"] > out["left_pos"])
        )
        cond_left_minus = (
            is_concordant
            & left_present
            & out["is_left"]
            & (ori0 == "-")
            & (out["query_pos"] < out["left_pos"])
        )
        drop_mask |= cond_left_plus | cond_left_minus

        right_present = out["right_pos"].notna()
        cond_right_plus = (
            is_concordant
            & right_present
            & ~(out["is_left"])
            & (ori1 == "+")
            & (out["query_end"] > out["right_pos"])
        )
        cond_right_minus = (
            is_concordant
            & right_present
            & ~(out["is_left"])
            & (ori1 == "-")
            & (out["query_pos"] < out["right_pos"])
        )
        drop_mask |= cond_right_plus | cond_right_minus
    else:
        out = output

    final_output = out.copy()
    if output_copy.shape[0] != 0:
        final_output["prune"] = drop_mask

    os.makedirs("fastq", exist_ok=True)
    for sv, reads in final_output.groupby(
        ["break_chrom1", "break_pos1", "break_chrom2", "break_pos2"]
    ):
        print("Making fastq for:", chrom1, pos1, chrom2, pos2)
        chrom1, pos1, chrom2, pos2 = sv
        fast1_name = f"fastq/b_{chrom1}_{pos1}_{chrom2}_{pos2}_1.fastq"
        fast2_name = f"fastq/b_{chrom1}_{pos1}_{chrom2}_{pos2}_2.fastq"
        with open(fast1_name, "w") as fast1:
            with open(fast2_name, "w") as fast2:
                it = reads.iterrows()
                for idx, read in it:
                    if "H" in read["query_cigar"]:
                        continue
                    if not output_copy.shape[0] == 0 and (read["read_num"] == "1") and (
                        (read["prune"]) or (reads.loc[idx + 1, "prune"])
                    ):
                        next(it)
                        continue
                    if read["read_num"] == "1":
                        fast1.write(
                            "@"
                            + read["query_name"]
                            + "\n"
                            + (
                                read["query_aln_full"]
                                if read["query_orientation"] == "+"
                                else rev_comp(read["query_aln_full"])
                            )
                            + "\n"
                            + "+\n"
                            + "".join(
                                list(
                                    map(
                                        lambda x: chr(x + 33),
                                        (
                                            read["query_qualities"]
                                            if read["query_orientation"] == "+"
                                            else read["query_qualities"][::-1]
                                        ),
                                    )
                                )
                            )
                            + "\n"
                        )
                    else:
                        fast2.write(
                            "@"
                            + read["query_name"]
                            + "\n"
                            + (
                                read["query_aln_full"]
                                if read["query_orientation"] == "+"
                                else rev_comp(read["query_aln_full"])
                            )
                            + "\n"
                            + "+\n"
                            + "".join(
                                list(
                                    map(
                                        lambda x: chr(x + 33),
                                        (
                                            read["query_qualities"]
                                            if read["query_orientation"] == "+"
                                            else read["query_qualities"][::-1]
                                        ),
                                    )
                                )
                            )
                            + "\n"
                        )
        subprocess.run(["gzip", fast1_name, "-f"])
        subprocess.run(["gzip", fast2_name, "-f"])

    final_output = final_output.drop(
        columns=["all_left_within_10", "all_right_within_10", "query_qualities"], errors="ignore"
    )
    leftover_output = leftover_output.drop(
        columns=["query_qualities"], errors="ignore"
    )

    if not has_homology:
        final_output = final_output.drop(columns=["homology_len", "homology_seq"], errors="ignore")
        leftover_output = leftover_output.drop(columns=["homology_len", "homology_seq"], errors="ignore")

    final_output.to_csv(
        args.file if args.file else args.bam.split("/")[-1].split(".")[0] + ".tsv",
        sep="\t",
        index=False,
    )
    leftover_output.to_csv(
        args.file if args.file else args.bam.split("/")[-1].split(".")[0] + "_leftover.tsv",
        sep="\t",
        index=False,
    )

    if args.verbose:
        print("Total Number of split reads:", num_split)
        print("Total Number of paired reads:", num_paired)
