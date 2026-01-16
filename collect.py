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


# @profile
# def get_pairs(
#     regions, chrom1, chrom2, break_pos1, break_pos2, sv_type, sv_orientation, args
# ):
#     """
#     Fetch all paired reads, including those with split alignments

#     Parameters:
#     - regions: the two SV regions to collect reads from
#     - chrom1: chromosome # of first region
#     - chrom2: chromosome # of second region

#     Returns:
#     - split_alignments: contains all sets of paired reads with split alignments
#     - paired_alignments: contains all sets of paired reads without split alignments
#     """
#     split_alignments = []
#     nonsplit_alignments = []

#     # Iterate over first breakpoint end
#     for i, region in enumerate(regions):
#         for aln in region:
#             # If we've seen this read before, skip it
#             if (
#                 aln.is_duplicate
#                 or any(aln in tup for tup in split_alignments)
#                 or any(aln in tup for tup in nonsplit_alignments)
#             ):
#                 # print("Duplicate read or read already seen")
#                 continue

#             if not (aln.is_mapped):
#                 print("Unmapped read found")
#                 # print(aln.query_sequence)
#                 # mate = samfile.mate(aln)
#                 # print(mate.reference_name)
#                 # print(mate.reference_start)
#                 # print("+" if mate.is_forward else "-")
#                 # print("read1" if mate.is_read1 else "read2")
#                 # print(mate.query_sequence)
#                 # print()
#                 continue

#             try:
#                 mate = samfile.mate(aln)
#             except ValueError as e:
#                 print(e)
#                 # print(aln.reference_name)
#                 # print(aln.reference_start)
#                 # print("+" if aln.is_forward else "-")
#                 # print("read1" if aln.is_read1 else "read2")
#                 # print(aln.query_sequence)
#                 # print()
#                 continue

#             # Filter out unhelpful discordant pairs
#             if not aln.is_proper_pair and not (
#                 aln.is_supplementary or aln.is_secondary or aln.has_tag("SA")
#             ):
#                 if (
#                     i == 0
#                     and (
#                         not mate.reference_name == chrom2
#                         or not (
#                             break_pos2 - args.refine < mate.reference_start
#                             and break_pos2 + args.refine > mate.reference_start
#                         )
#                     )
#                 ) or (
#                     i == 1
#                     and (
#                         not mate.reference_name == chrom1
#                         or not (
#                             break_pos1 - args.refine < mate.reference_start
#                             and break_pos1 + args.refine > mate.reference_start
#                         )
#                     )
#                 ):
#                     continue

#             # Avoid inclusion of non-informative region to the right of foldback
#             # if sv_type == "foldback":
#             #     if (
#             #         i == 0
#             #         and (
#             #             aln.reference_end <= break_pos1 - 20
#             #             or (
#             #                 mate.reference_name == chrom1
#             #                 and mate.reference_end <= break_pos1 - 20
#             #             )
#             #         )
#             #         or (
#             #             mate.reference_name == chrom2
#             #             and mate.reference_start > break_pos2 + 20
#             #         )
#             #     ):
#             #         continue
#             #     if (
#             #         i == 1
#             #         and (
#             #             aln.reference_start > break_pos2 + 20
#             #             or (
#             #                 mate.reference_name == chrom2
#             #                 and mate.reference_start > break_pos2 + 20
#             #             )
#             #         )
#             #         or (
#             #             mate.reference_name == chrom1
#             #             and mate.reference_end <= break_pos1 - 20
#             #         )
#             #     ):
#             #         continue
#             read1, sup, read2 = None, None, None
#             if aln.is_supplementary or aln.is_secondary or aln.has_tag("SA"):
#                 sa_tag = (aln.get_tag("SA")).split(";")[:-1]
#                 sa_tag = sa_tag[0]
#                 sup_name, sup_pos, sup_strand, sup_cigar, sup_mapq, sup_tlen = (
#                     sa_tag.split(",")
#                 )
#                 sup_pos = int(sup_pos)
#                 sup_tlen = int(sup_tlen)

#                 found_split2 = False
#                 for sup_read in samfile.fetch(sup_name, sup_pos, sup_pos + 1):
#                     if (
#                         sup_read.query_name == aln.query_name
#                         and sup_read.reference_name == sup_name
#                         and (sup_read.reference_start + 1) == sup_pos
#                         and _sh_re.sub("", sup_read.cigarstring)
#                         == _sh_re.sub("", sup_cigar)
#                     ):
#                         non_sup = (
#                             aln
#                             if not aln.is_supplementary and not aln.is_secondary
#                             else sup_read
#                         )
#                         sup = (
#                             aln
#                             if aln.is_supplementary or aln.is_secondary
#                             else sup_read
#                         )
#                         read1 = non_sup if non_sup.is_read1 else mate
#                         read2 = non_sup if non_sup.is_read2 else mate
#                         found_split2 = True
#                         break
#                 if not found_split2:
#                     print(sup_name, sup_pos, sup_cigar, sup_tlen)
#                     raise NameError("split mate not found")

#                 if not filter_for_pos_ori(
#                     read1,
#                     read2,
#                     chrom1,
#                     break_pos1,
#                     chrom2,
#                     break_pos2,
#                     sv_type,
#                     sv_orientation,
#                     sup,
#                 ):
#                     continue

#                 split_alignments.append((read1, sup, read2))
#                 split_alignments.append((read1, sup, read2))

#             elif mate.is_supplementary or mate.is_secondary or mate.has_tag("SA"):
#                 sa_tag = (mate.get_tag("SA")).split(";")[:-1]
#                 sa_tag = sa_tag[0]
#                 sup_name, sup_pos, sup_strand, sup_cigar, sup_mapq, sup_tlen = (
#                     sa_tag.split(",")
#                 )
#                 sup_pos = int(sup_pos)
#                 sup_tlen = int(sup_tlen)

#                 found_split2 = False
#                 for sup_read in samfile.fetch(sup_name, sup_pos, sup_pos + 1):
#                     if (
#                         sup_read.query_name == mate.query_name
#                         and sup_read.reference_name == sup_name
#                         and (sup_read.reference_start + 1) == sup_pos
#                         and _sh_re.sub("", sup_read.cigarstring)
#                         == _sh_re.sub("", sup_cigar)
#                     ):
#                         non_sup = (
#                             mate
#                             if not mate.is_supplementary and not mate.is_secondary
#                             else sup_read
#                         )
#                         sup = (
#                             mate
#                             if mate.is_supplementary or mate.is_secondary
#                             else sup_read
#                         )
#                         read1 = non_sup if non_sup.is_read1 else aln
#                         read2 = non_sup if non_sup.is_read2 else aln
#                         found_split2 = True
#                         break

#                 if not found_split2:
#                     print(sup_name, sup_pos, sup_cigar, sup_tlen)
#                     raise NameError("split mate not found")

#                 if not filter_for_pos_ori(
#                     read1,
#                     read2,
#                     chrom1,
#                     break_pos1,
#                     chrom2,
#                     break_pos2,
#                     sv_type,
#                     sv_orientation,
#                     sup,
#                 ):
#                     continue

#                 split_alignments.append((read1, sup, read2))
#                 split_alignments.append((read1, sup, read2))

#             elif not aln.is_proper_pair and aln.next_reference_name == (
#                 chrom2 if i == 0 else chrom1
#             ):
#                 read1 = aln if aln.is_read1 else mate
#                 read2 = aln if aln.is_read2 else mate

#                 if not filter_for_pos_ori(
#                     read1,
#                     read2,
#                     chrom1,
#                     break_pos1,
#                     chrom2,
#                     break_pos2,
#                     sv_type,
#                     sv_orientation,
#                 ):
#                     continue

#                 nonsplit_alignments.append((read1, read2))
#                 nonsplit_alignments.append((read1, read2))

#             elif aln.is_proper_pair:
#                 read1 = aln if aln.is_read1 else mate
#                 read2 = aln if aln.is_read2 else mate
#                 nonsplit_alignments.append((read1, read2))

#     return split_alignments, nonsplit_alignments


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

    # print(reads_by_name)
    # print("HWI-ST1113:353:H99R3ADXX:2:1106:18680:72982" in reads_by_name)
    # print(len(reads_by_name["HWI-ST1113:353:H99R3ADXX:2:1106:18680:72982"]))
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
                    else:
                        split_leftover_alignments.append((read1, sup, read2))

                else:
                    # print(qname)
                    # print(sup_or_sec[0].reference_name, sup_or_sec[1].reference_name)
                    # print(sup_or_sec[0].reference_start, sup_or_sec[1].reference_start)
                    # print()

                    sup1, sup2 = (
                        sup_or_sec[0],
                        sup_or_sec[1]
                    ) if sup_or_sec[0].is_read1 else (sup_or_sec[1], sup_or_sec[0])

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
                        split_leftover_alignments.append((read1, sup1, read2, sup2))

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
    amplicon,
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
    # print(chrom1, pos1, args.refine)
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
                            "amplicon": amplicon,
                            "query_qualities": read.query_qualities,
                        }
                    )
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
                            "amplicon": amplicon,
                            "query_qualities": read.query_qualities,
                        }
                    )
    
    split_leftover_alignment_details = []
    for pair in split_leftover_alignments:
        if len(pair) == 3:
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
                    query_aln_full = read.query_sequence
                    if read.is_secondary or read.is_supplementary:
                        prim_read = read1_1 if read1_1.has_tag("SA") else read2
                        query_aln_full = (
                            prim_read.query_sequence
                            if prim_read.is_forward == read.is_forward
                            else rev_comp(prim_read.query_sequence)
                        )

                    split_leftover_alignment_details.append(
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
                            "amplicon": amplicon,
                            "query_qualities": read.query_qualities,
                        }
                    )
        if len(pair) == 4:
            read1_1, read1_2, read2_1, read2_2 = pair
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

                    split_leftover_alignment_details.append(
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
                            "amplicon": amplicon,
                            "query_qualities": read.query_qualities,
                        }
                    )

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
                        "amplicon": amplicon,
                        "query_qualities": read.query_qualities,
                    }
                )

    # Convert details to pandas DataFrames
    split_df = pd.DataFrame(split_alignment_details)
    nonsplit_df = pd.DataFrame(nonsplit_alignment_details)
    split_leftover_df = pd.DataFrame(split_leftover_alignment_details)
    # paired_df = pd.DataFrame(paired_alignment_details)

    return split_df, nonsplit_df, split_leftover_df


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
                amplicon=f.split("_")[1]
            )
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

    output = []
    leftover_output = []
    for index, row in enumerate(df.itertuples(index=False), 1):
        print(f"Breakpoint {index}")
        chrom1 = row.chrom1
        pos1 = row.pos1
        chrom2 = row.chrom2
        pos2 = row.pos2
        sv_type = row.sv_type
        read_support = row.read_support
        features = row.features
        orientation = row.orientation
        hom_len = row.homology_length
        homology = row.homology_sequence
        amplicon = row.amplicon
        # print(chrom1, pos1, chrom2, pos2)

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
            amplicon,
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
        # print(drop_mask.iloc[406])

        right_present = out["right_pos"].notna()
        cond_right_plus = (
            is_concordant
            & right_present
            & ~(out["is_left"])
            & (ori1 == "+")
            & (out["query_end"] > out["right_pos"])
        )
        # print(cond_right_plus.iloc[406])
        cond_right_minus = (
            is_concordant
            & right_present
            & ~(out["is_left"])
            & (ori1 == "-")
            & (out["query_pos"] < out["right_pos"])
        )
        # print(cond_right_minus.iloc[406])
        # print(out["query_pos"].iloc[406], out["right_pos"].iloc[406])
        drop_mask |= cond_right_plus | cond_right_minus
        # print(drop_mask.iloc[406])
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
