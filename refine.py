#!/usr/bin/env python3
import os
import argparse
import warnings
import pandas as pd
from Bio import Align
import subprocess, shutil
from natsort import index_natsorted
import numpy as np
import marisa_trie
import math
import re

# suppress pandas warnings
warnings.filterwarnings("ignore")
pd.set_option("mode.chained_assignment", None)

bp_to_read_idxs = {}

# -------------------------------------------
# Utility functions
# -------------------------------------------

trans = str.maketrans("ACGT", "TGCA")

def rev_comp_vec(arr):
    return np.frompyfunc(lambda s: s[::-1].translate(trans), 1, 1)(arr).astype(str)

def find_longest_common_variant(strs):
    trie = marisa_trie.Trie(strs.tolist())
    result = None
    count = 0
    for strn in strs:
        num_prefixes = len(trie.prefixes(strn))
        if num_prefixes > count:
            count = num_prefixes
            result = strn
    return result

def top3_find_longest_common_variant(strs):
    strs_list = strs.tolist()
    str_counts = dict({ k:0 for k in strs_list })
    trie = marisa_trie.Trie(strs_list)
    for strn in strs:
        num_prefixes = len(trie.prefixes(strn))
        str_counts[strn] += num_prefixes if strn not in str_counts else 1
    sorted_str_counts = list(sorted(str_counts.items(), key = lambda item: item[1]))
    return np.array([seq for seq, val in sorted_str_counts[-3:]])

def rev_comp(seq):
    trans = str.maketrans("ACGT", "TGCA")
    return (
        seq[::-1].translate(trans)
        if isinstance(seq, str)
        else seq.str[::-1].str.translate(trans)
    )

def get_homology(first, last):
    for i in range(len(first), -1, -1):
        if last.startswith(first[-i:]):
            return first[-i:]
    return ""


# -------------------------------------------
# Refine pipeline
# -------------------------------------------

print_cols = [
    "query_name",
    "proper_pair",
    "query_pos",
    "query_end",
    "query_cigar",
    "split",
    "query_aln_full",
]
print_cols2 = [
    "query_name",
    "query_short",
    "query_chrom",
    "query_pos",
    "query_cigar",
    "hom_clip_match",
    "ins_clip_match",
    "split",
    "rev_clipped",
    "query_aln_sub",
]

def refine_step1(reads, all_reads, verbose, is_left):
    reads = reads[reads["query_cigar"].str.contains(r"[SH]")].copy()
    reads["begin"] = reads["query_cigar"].str.contains(r"^\d+[SH]", regex=True)
    leftover["begin"] = leftover["query_cigar"].str.contains(r"^\d+[SH]", regex=True)
    reads["end"] = reads["query_cigar"].str.contains(r"\d+[SH]$", regex=True)
    leftover["end"] = leftover["query_cigar"].str.contains(r"\d+[SH]$", regex=True)
    reads["sv_end"] = reads["query_pos"].where(reads["begin"])
    reads["sv_end"].fillna(reads["query_end"].where(reads["end"]), inplace=True)
    reads.dropna(subset=["sv_end"], inplace=True)
    mask_ori = ((reads["break_orientation"].str.get(0) == "+") & reads["end"]) | ((reads["break_orientation"].str.get(0) == "-") & reads["begin"]) if is_left else ((reads["break_orientation"].str.get(1) == "+") & reads["end"]) | ((reads["break_orientation"].str.get(1) == "-") & reads["begin"])
    mask_not_mul_clips = (reads["begin"] ^ reads["end"])
    reads = reads[(mask_ori & mask_not_mul_clips)]

    if verbose >= 2:
        print(
            reads.sort_values("sv_end")[["sv_end", *print_cols]].to_string(index=False),
            "\n",
        )

    cands = reads["sv_end"].value_counts()
    best = cands.index[0] if len(cands) > 0 else None

    if verbose == 1 and best is not None:
        sel = reads[reads["sv_end"] == best][["sv_end", *print_cols]]
        disc = sel[sel["proper_pair"] == "Discordant"]
        mates = all_reads[
            all_reads["query_name"].isin(disc.query_name)
            & ~all_reads["query_aln_full"].isin(disc.query_aln_full)
        ][print_cols]
        print(
            sel.to_string(index=False), "\nMates:\n", mates.to_string(index=False), "\n"
        )

    reads['query_name'] = np.where(reads['read_num']==1,
                      reads['query_name'] + '_1',
                      np.where(reads['read_num']==2,
                               reads['query_name'] + '_2',
                               reads['query_name']))

    return best, reads

def in_range(test, middle, width):
    if test >= middle - width and test <= middle + width:
        return True
    return False
def cord_sort(row1, row2):
    def map_to_int(s):
        return 22 if s == "X" else 23 if s == "Y" else int(s)
    row1_chr, row1_pos = row1[["query_chrom", "query_pos"]]
    row1_chr = map_to_int(row1_chr[3:])
    row2_chr, row2_pos = row2[["query_chrom", "query_pos"]]
    row2_chr = map_to_int(row2_chr[3:])
    rowL, rowR = (row1, row2)
    if row1_chr > row2_chr:
        rowL, rowR = (row2, row1)
    elif row1_chr == row2_chr:
        if row1_pos > row2_pos:
            rowL, rowR = (row2, row1)
    
    return rowL, rowR

def find_top_variant(left, right):
    fr = find_longest_common_variant(left["query_aln_sub"]) if left["break_orientation"].iloc[0][0] == "-" else rev_comp(find_longest_common_variant(rev_comp_vec(left["query_aln_sub"])))
    mask = left["query_aln_sub"].to_numpy() == fr
    idx = np.flatnonzero(mask)[0]
    fr = left.iloc[idx]

    lr = find_longest_common_variant(right["query_aln_sub"]) if right["break_orientation"].iloc[0][1] == "-" else rev_comp(find_longest_common_variant(rev_comp_vec(right["query_aln_sub"])))
    mask = right["query_aln_sub"].to_numpy() == lr
    idx = np.flatnonzero(mask)[0]
    lr = right.iloc[idx]

    l_ori = fr.break_orientation[0] if fr.begin and fr.end else ("+" if fr.end else "-")
    r_ori = lr.break_orientation[1] if lr.begin and lr.end else ("+" if lr.end else "-")
    
    return fr, lr, l_ori, r_ori

def find_top_insertion(left, right, seq1=None, seq2=None, ins=None):
    c1 = find_longest_common_variant(left["rev_clipped"]) if left["break_orientation"].iloc[0][0] == "-" else find_longest_common_variant(left["clipped"])
    mask = left["rev_clipped"].to_numpy() == c1
    idx = np.flatnonzero(mask)[0]
    c1 = left.iloc[idx]["rev_clipped"]

    c2 = rev_comp(find_longest_common_variant(rev_comp_vec(right["clipped"]))) if right["break_orientation"].iloc[0][1] == "-" else rev_comp(find_longest_common_variant(right["clipped"]))
    mask = right["rev_clipped"].to_numpy() == c2
    idx = np.flatnonzero(mask)[0]
    c2 = right.iloc[idx]["rev_clipped"]
    if not ins:
        ins = get_homology(c2, c1)
    ilen = len(ins)

    if not seq1 and not seq2:
        fr, lr, l_ori, r_ori = find_top_variant(left, right)
        seq1 = fr.query_aln_sub if l_ori == "+" else rev_comp(fr.query_aln_sub)
        seq2 = lr.query_aln_sub if r_ori == "-" else rev_comp(lr.query_aln_sub)
    
    def check_starts(seq):
        return seq2.startswith(seq[ilen:]) if len(seq) != 0 else True
    def check_ends(seq):
        return (seq1.endswith(seq[:-ilen]) if len(seq) != 0 else True) if ilen != 0 else False
    
    left["ins_clip_match"] = np.frompyfunc(check_starts, 1, 1)(left["rev_clipped"])
    right["ins_clip_match"] = np.frompyfunc(check_ends, 1, 1)(right["rev_clipped"])

    isum_l = left["ins_clip_match"].sum()
    isum_r = right["ins_clip_match"].sum()

    return ins, ilen, isum_l, isum_r

def get_templated_ins(read_first, read_last, tst_reads):
    template = read_first["query_aln_full"]
    read1_offset = int(re.search(r"^\d+[M]", read_first["query_cigar"]).group(0)[:-1])
    read2_offset = int(re.search(r"^\d+[SH]", read_last["query_cigar"]).group(0)[:-1]) if read_last["query_orientation"]==read_first["query_orientation"] else int(re.search(r"\d+[SH]$", read_last["query_cigar"]).group(0)[:-1])
    
    tst = template[read1_offset:read2_offset]
    ref_cords = []
    prev_offset = read1_offset
    ref_start = None
    ref_end = None
    is_same_ori = False
    for idx, read in tst_reads.iterrows():
        is_same_ori = read["query_orientation"]==read_first["query_orientation"]
        read_offset = int(re.search(r"^\d+[SH]", read["query_cigar"]).group(0)[:-1]) if is_same_ori else int(re.search(r"\d+[SH]$", read["query_cigar"]).group(0)[:-1])
        if idx != 0:
            if is_same_ori:
                ref_cords.append(ref_end - (prev_offset - read_offset))
            else:
                ref_cords.append(ref_start + (prev_offset - read_offset))
        ref_start, ref_end = read[["query_pos", "query_end"]]
        ref_cords.append(read["query_chrom"])
        if is_same_ori:
            ref_cords.append(ref_start + (prev_offset - read_offset))
        else:
            ref_cords.append(ref_end - (prev_offset - read_offset))

        ref_cords.append(read["query_chrom"])
        prev_offset += len(read["query_aln_sub"]) - (prev_offset - read_offset)
    if is_same_ori:
        ref_cords.append(ref_end - (prev_offset - read2_offset))
    else:
        ref_cords.append(ref_start + (prev_offset - read2_offset))

    return tst, ref_cords

# @profile
def check_overlap(left, right, leftover):
    results = []

    # compute clipped sequences
    for df in (left, right, leftover):
        df["clip_len"] = df["query_aln_full"].str.len() - df["query_aln_sub"].str.len()
        if not df.empty:
            df["clipped"] = df.apply(
                lambda x: (
                    x["query_aln_full"][: x["clip_len"]]
                    if x["begin"]
                    else x["query_aln_full"][
                        (
                            -x["clip_len"]
                            if x["clip_len"] <= len(x["query_aln_full"])
                            else len(x["query_aln_full"])
                        ) :
                    ]
                ),
                axis=1,
            )
        else:
            df["clipped"] = pd.Series(dtype=str)

    # # pick the longest sub-alignment from each
    lmask = (~left["end"]).to_numpy()
    rmask = (~right["begin"]).to_numpy()

    left["rev_clipped"] = np.where(lmask, rev_comp_vec(left["clipped"].to_numpy().astype(str)), left["clipped"])
    right["rev_clipped"] = np.where(rmask, rev_comp_vec(right["clipped"].to_numpy().astype(str)), right["clipped"])

    cand_svs = {}
    break_pos1 = left.iloc[0]["break_pos1"]
    if break_pos1 in bp_to_read_idxs:
        for tst_idxs in bp_to_read_idxs[break_pos1]:

            if len(tst_idxs) < 3:
                continue

            mask = (left["break_pos1"] == break_pos1)
            idx = np.flatnonzero(mask)[0]
            break_pos2, break_chr1, break_chr2, break_ori = left.iloc[idx][["break_pos2", "break_chrom1", "break_chrom2", "break_orientation"]]
            read_first = leftover.loc[tst_idxs[0]]
            read_last = leftover.loc[tst_idxs[-1]]
            read_l, read_r = cord_sort(read_first, read_last)
            ori_l = "query_pos" if break_ori[0] == "-" else "query_end"
            lv = read_l[ori_l]
            chr_l = read_l["query_chrom"]
            ori_r = "query_pos" if break_ori[1] == "-" else "query_end"
            rv = read_r[ori_r]
            chr_r = read_r["query_chrom"]
            if chr_l != break_chr1 or chr_r != break_chr2 or not in_range(lv, break_pos1, 500) or not in_range(rv, break_pos2, 500):
                continue
            if (lv, rv, chr_l, chr_r) in cand_svs or (break_ori[0] == "-" and read_l["end"]) or (break_ori[0] == "+" and read_l["begin"]) or (break_ori[1] == "-" and read_r["end"]) or (break_ori[1] == "+" and read_r["begin"]):
                continue
            
            lmask = (leftover["query_chrom"] == chr_l) & (leftover[ori_l] == lv)
            rmask = (leftover["query_chrom"] == chr_r) & (leftover[ori_r] == rv)
            leftover_left = leftover[lmask].copy()
            leftover_right = leftover[rmask].copy()
            lmask = (~leftover_left["end"]).to_numpy()
            rmask = (~leftover_right["begin"]).to_numpy()

            leftover_left["rev_clipped"] = np.where(lmask, rev_comp_vec(leftover_left["clipped"].to_numpy().astype(str)), leftover_left["clipped"])
            leftover_right["rev_clipped"] = np.where(rmask, rev_comp_vec(leftover_right["clipped"].to_numpy().astype(str)), leftover_right["clipped"])
            tst, ref_cords = get_templated_ins(read_first, read_last, leftover.loc[tst_idxs[1:-1]].reset_index(drop=True))
            cand_svs[(lv, rv, chr_l, chr_r)] = [leftover_left, leftover_right, tst, ref_cords]

        for cand_sv, value in cand_svs.items():
            lv, rv, chr_l, chr_r = cand_sv
            leftover_left, leftover_right, tst, ref_cords = value

            left_copy = left[
                ~left["split"]
                | (
                    left["split"]
                    & left["query_name"].isin(right.loc[right["split"], "query_name"])
                )
            ]
            right_copy = right[
                ~right["split"]
                | (
                    right["split"]
                    & right["query_name"].isin(left_copy.loc[left_copy["split"], "query_name"])
                )
            ]

            lmask = (left_copy["query_chrom"] == chr_l) & (left_copy["query_pos" if break_ori[0] == "-" else "query_end"] == lv)
            rmask = (right_copy["query_chrom"] == chr_r) & (right_copy["query_pos" if break_ori[1] == "-" else "query_end"] == rv)
            left_copy = pd.concat([left_copy[lmask], leftover_left], join='inner')
            right_copy = pd.concat([right_copy[rmask], leftover_right], join='inner')
            left_copy["hom_clip_match"] = False
            right_copy["hom_clip_match"] = False

            _, _, isum_l, isum_r = find_top_insertion(left_copy, right_copy, ins=tst)
            fr, lr, l_ori, r_ori = find_top_variant(left_copy, right_copy)

            total_reads = len(left_copy) + len(right_copy)
            num_split_reads = left_copy["split"].sum()

            lmask = (~left_copy["split"])
            rmask = (~right_copy["split"])

            left_matches_ins = left_copy["ins_clip_match"].where(lmask).sum()
            right_matches_ins = right_copy["ins_clip_match"].where(rmask).sum()
            
            num_left_soft_reads = (~left_copy["split"]).sum()
            num_right_soft_reads = (~right_copy["split"]).sum()
            results.append(
                pd.DataFrame(
                    {
                        "left_sv": [lv],
                        "right_sv": [rv],
                        "hom_sum_left": [0],
                        "hom_sum_right": [0],
                        "ins_sum_left": [isum_l],
                        "ins_sum_right": [isum_r],
                        "homology": [""],
                        "insertion": [tst],
                        "tst_cords": [ref_cords],
                        "left_len": [len(left_copy)],
                        "right_len": [len(right_copy)],
                        "split_support": [num_split_reads],
                        "soft_support": [num_left_soft_reads + num_right_soft_reads],
                        "left_soft_matches_hom": [0],
                        "right_soft_matches_hom": [0],
                        "left_soft_matches_ins": [left_matches_ins],
                        "right_soft_matches_ins": [right_matches_ins],
                        "total_reads": [total_reads],
                        "lcommon_var": [fr["query_aln_sub"]],
                        "rcommon_var": [lr["query_aln_sub"]],
                        # <-- store DataFrame directly, not as tuple
                        "left": [left_copy.copy()],
                        "right": [right_copy.copy()],
                    }
                )
            )

    for lv, ldf in left.groupby("sv_end"):
        for rv, rdf in right.groupby("sv_end"):

            # filter splits
            ldf_copy = ldf.copy()
            ldf_copy = ldf_copy[
                ~ldf_copy["split"]
                | (
                    ldf_copy["split"]
                    & ldf_copy["query_name"].isin(rdf.loc[rdf["split"], "query_name"])
                )
            ]
            rdf = rdf[
                ~rdf["split"]
                | (
                    rdf["split"]
                    & rdf["query_name"].isin(ldf_copy.loc[ldf_copy["split"], "query_name"])
                )
            ]

            if ldf_copy.empty or rdf.empty or len(ldf_copy) < 3 or len(rdf) < 3:
                continue

            fr = top3_find_longest_common_variant(ldf_copy["query_aln_sub"]) if ldf_copy["break_orientation"].iloc[0][0] == "-" else [rev_comp(x) for x in top3_find_longest_common_variant(ldf_copy["query_aln_sub"].map(rev_comp))]
            mask = np.isin(ldf_copy["query_aln_sub"].to_numpy(), fr)
            idx = np.flatnonzero(mask)
            fr = ldf_copy.iloc[idx].drop_duplicates(subset="query_aln_sub")

            lr = top3_find_longest_common_variant(rdf["query_aln_sub"]) if rdf["break_orientation"].iloc[0][1] == "-" else [rev_comp(x) for x in top3_find_longest_common_variant(rdf["query_aln_sub"].map(rev_comp))]
            mask = np.isin(rdf["query_aln_sub"].to_numpy(), lr)
            idx = np.flatnonzero(mask)
            lr = rdf.iloc[idx].drop_duplicates(subset="query_aln_sub")

            for row_l in fr.itertuples():
                for row_r in lr.itertuples():
                    lcommon_var = row_l.rev_clipped
                    rcommon_var = row_r.rev_clipped
                    l_ori = row_l.break_orientation[0] if row_l.begin and row_l.end else ("+" if row_l.end else "-")
                    r_ori = row_r.break_orientation[1] if row_r.begin and row_r.end else ("+" if row_r.end else "-")

                    seq1 = row_l.query_aln_sub if l_ori == "+" else rev_comp(row_l.query_aln_sub)
                    seq2 = row_r.query_aln_sub if r_ori == "-" else rev_comp(row_r.query_aln_sub)

                    hom = get_homology(seq1, seq2)
                    hlen = len(hom)
                    s1_no = seq1[:-hlen] if hlen != 0 else seq1
                    s2_no = seq2[hlen:]

                    # homology matching
                    for df, target_no, name in (
                        (ldf_copy, s2_no, "hom_clip_match"),
                        (rdf, s1_no, "hom_clip_match"),
                    ):
                        def check_starts(seq):
                            return target_no.startswith(seq) or target_no in seq
                        def check_ends(seq):
                            return target_no.endswith(seq) or target_no in seq
                        if df is ldf_copy:
                            df[name] = np.frompyfunc(check_starts, 1, 1)(df["rev_clipped"])
                        else:
                            df[name] = np.frompyfunc(check_ends, 1, 1)(df["rev_clipped"])
                    
                    hsum_l = ldf_copy["hom_clip_match"].sum()
                    hsum_r = rdf["hom_clip_match"].sum()

                    # insertion matching
                    ins, ilen, isum_l, isum_r = find_top_insertion(ldf_copy, rdf, seq1, seq2)

                    if hsum_l == 0 and hsum_r == 0 and hlen != 0:
                        for df, target_no, name in (
                        (ldf_copy, seq2, "hom_clip_match"),
                        (rdf, seq1, "hom_clip_match"),
                        ):

                            def check_starts(seq):
                                return target_no.startswith(seq) or target_no in seq
                            def check_ends(seq):
                                return target_no.endswith(seq) or target_no in seq
                            if df is ldf:
                                df[name] = np.frompyfunc(check_starts, 1, 1)(df["rev_clipped"])
                            else:
                                df[name] = np.frompyfunc(check_ends, 1, 1)(df["rev_clipped"])

                        hsum_l = ldf_copy["hom_clip_match"].sum()
                        hsum_r = rdf["hom_clip_match"].sum()
                        hom = ""
                        hlen = 0
                
                    
                    if (hsum_l == 0 or hsum_r == 0) and (isum_l == 0 or isum_r == 0):
                        continue

                    total_reads = len(ldf_copy) + len(rdf)
                    num_split_reads = ldf_copy["split"].sum()

                    lmask = (~ldf_copy["split"])
                    rmask = (~rdf["split"])

                    left_matches_hom = ldf_copy["hom_clip_match"].where(lmask).sum()
                    right_matches_hom = rdf["hom_clip_match"].where(rmask).sum()
                    left_matches_ins = ldf_copy["ins_clip_match"].where(lmask).sum()
                    right_matches_ins = rdf["ins_clip_match"].where(rmask).sum()
                    
                    num_left_soft_reads = (~ldf_copy["split"]).sum()
                    num_right_soft_reads = (~rdf["split"]).sum()

                    if num_split_reads == 0:
                        hom_left = hsum_l / num_left_soft_reads
                        hom_right = hsum_r / num_right_soft_reads
                        ins_left = isum_l / num_left_soft_reads
                        ins_right = isum_r / num_right_soft_reads

                        if not (((hom_left >= 0.33) and (hom_right >= 0.33)) or ((ins_left >= 0.33) and (ins_right >= 0.33))):
                            continue

                    results.append(
                        pd.DataFrame(
                            {
                                "left_sv": [lv],
                                "right_sv": [rv],
                                "hom_sum_left": [hsum_l],
                                "hom_sum_right": [hsum_r],
                                "ins_sum_left": [isum_l],
                                "ins_sum_right": [isum_r],
                                "homology": [hom],
                                "insertion": [ins],
                                "tst_cords": [],
                                "left_len": [len(ldf_copy)],
                                "right_len": [len(rdf)],
                                "split_support": [num_split_reads],
                                "soft_support": [num_left_soft_reads + num_right_soft_reads],
                                "left_soft_matches_hom": [left_matches_hom],
                                "right_soft_matches_hom": [right_matches_hom],
                                "left_soft_matches_ins": [left_matches_ins],
                                "right_soft_matches_ins": [right_matches_ins],
                                "total_reads": [total_reads],
                                "lcommon_var": [lcommon_var],
                                "rcommon_var": [rcommon_var],
                                "left": [ldf_copy.copy()],
                                "right": [rdf.copy()],
                            }
                        )
                    )
    if not results:
        return None
    out = pd.concat(results, ignore_index=True)
    out["hom_%"] = (out["hom_sum_left"] + out["hom_sum_right"]) / (
        out["left_len"] + out["right_len"]
    )
    out["ins_%"] = (out["ins_sum_left"] + out["ins_sum_right"]) / (
        out["left_len"] + out["right_len"]
    )
    out["total_%"] = out[["hom_%", "ins_%"]].max(axis=1)

    return out

def run_split(args):
    all_reads = pd.read_csv(args.file, sep="\t").drop_duplicates()
    leftover_splits = pd.read_csv(args.file[:-4] + "_leftover" + ".tsv", sep="\t").drop_duplicates()
    has_homology = ("AA_homology_len" in all_reads.columns) and ("AA_homology_seq" in all_reads.columns)

    for group, reads in leftover_splits.groupby(["break_pos1", "query_name", "read_num"]):
        if len(reads) > 2:
            bp_to_read_idxs.setdefault(group[0], []).append(reads.index.to_list())
    svs = all_reads.groupby("break_pos1")

    split_log = open(args.split_log + ".txt", "w")
    summary = []

    if args.list:
        cols = ["break_chrom1", "break_pos1", "break_chrom2", "break_pos2"]
        if has_homology:
            cols.append("AA_homology_seq")
        tbl = svs[cols].first().reset_index(drop=True)
        tbl.index.name = "breakpoint_index"
        print(tbl.to_string())
        return

    picks = None
    if args.breakpoints:
        idxs = svs["break_pos1"].first().reset_index(drop=True)
        picks = set(idxs.iloc[args.breakpoints].tolist())

    for bp, sv in svs:
        if picks and bp not in picks:
            continue

        if sv["break_chrom1"].iloc[0] == sv["break_chrom2"].iloc[0]:
            left = sv[
                (sv["query_pos"] >= bp - 500)
                & (sv["query_end"] <= bp + 500)
                & (sv["query_end"] < sv["break_pos2"].iloc[0] - 150)
            ]
            right = sv[
                (sv["query_pos"] >= sv["break_pos2"].iloc[0] - 500)
                & (sv["query_end"] <= sv["break_pos2"].iloc[0] + 500)
                & (sv["query_pos"] > bp + 150)
            ]
        else:
            left = sv[(sv["query_pos"] >= bp - 500) & (sv["query_end"] <= bp + 500)]
            right = sv[
                (sv["query_pos"] >= sv["break_pos2"].iloc[0] - 500)
                & (sv["query_end"] <= sv["break_pos2"].iloc[0] + 500)
            ]

        _, lgrp = refine_step1(left, all_reads, leftover_splits, args.verbose, True)
        _, rgrp = refine_step1(right, all_reads, leftover_splits, args.verbose, False)

        res = check_overlap(lgrp, rgrp, leftover_splits)

        if res is None or len(res) == 0:
            print("No overlap found")
            continue

        row = sv.iloc[0]
        break_pos1 = row["break_pos1"]
        break_pos2 = row["break_pos2"]
        break_ori = row["break_orientation"]
        left_within_20 = (res["left_sv"] >= break_pos1-20) & (res["left_sv"] <= break_pos1+20)
        right_within_20 = (res["right_sv"] >= break_pos2-20) & (res["right_sv"] <= break_pos2+20)
        left_match_ori = (res["left_sv"] <= break_pos1) if break_ori[0] == "-" else (res["left_sv"] >= break_pos1)
        right_match_ori = (res["right_sv"] <= break_pos2) if break_ori[1] == "-" else (res["right_sv"] >= break_pos2)
        res_mask = (left_within_20 | left_match_ori) & (right_within_20 | right_match_ori)
        res = res[res_mask]
        res = res.sort_values(["split_support", "total_%", "total_reads"], ascending=False).drop_duplicates(subset=["left_sv", "right_sv", "split_support", "total_%", "homology", "insertion"])

        res = res.head(3)

        if res is None or len(res) == 0:
            print("No overlap found")
            continue

        pd.set_option('display.max_colwidth', None)
        split_log.write(f"=== Breakpoint {bp} ===\n")
        for idx, row in res.iterrows():
            split_log.write(f"-- Candidate {idx} --\n")
            split_log.write(row.drop(["left", "right"]).to_string() + "\n")
            split_log.write(
                "Left reads:\n" + row["left"][print_cols2].to_string() + "\n"
            )
            split_log.write(
                "Right reads:\n" + row["right"][print_cols2].to_string() + "\n\n"
            )

        if args.breakpoints:
            for i in range(len(res) - 1, -1, -1):
                row = res.iloc[i]
                print("_" * 80)
                print(row.drop(["left", "right"]).to_string())
                print("Left reads:\n", row["left"][print_cols2].to_string())
                print("Right reads:\n", row["right"][print_cols2].to_string())
                print("_" * 80, "\n")

        top = res.iloc[0]
        summary.append(
            {
                "break_pos1": bp,
                "split_support": top.split_support,
                "soft_support": top.soft_support,
                "left_soft_matches": top.left_soft_matches_hom
                if top["hom_%"] > top["ins_%"]
                else top.left_soft_matches_ins,
                "right_soft_matches": top.right_soft_matches_hom
                if top["hom_%"] > top["ins_%"]
                else top.right_soft_matches_ins,
                "sp_left_sv": int(top.left_sv),
                "sp_right_sv": int(top.right_sv),
                "sp_hom_len": len(top.homology)
                if top["hom_%"] > top["ins_%"]
                else -len(top.insertion) if isinstance(top.insertion, str) else "N/A",
                "hom": top.homology if top["hom_%"] > top["ins_%"] else top.insertion,
            }
        )

        if args.breakpoints and has_homology:
            orig = sv.iloc[0]
            print(
                "Original AA breakpoint:\n",
                orig[
                    [
                        "break_chrom1",
                        "break_pos1",
                        "break_chrom2",
                        "break_pos2",
                        "AA_homology_seq",
                        "break_orientation",
                    ]
                ].to_string(),
                "\n",
            )

    split_log.close()

    brk_cols = [
        "break_chrom1",
        "break_pos1",
        "break_chrom2",
        "break_pos2",
        "break_sv_type",
        "break_read_support",
        "break_features",
        "break_orientation",
        "amplicon",
    ]
    if has_homology:
        brk_cols.insert(8, "AA_homology_len")
        brk_cols.insert(9, "AA_homology_seq")

    brk = svs[brk_cols].first().reset_index(drop=True)

    aug = pd.DataFrame(summary, columns=["break_pos1","split_support","soft_support","left_soft_matches","right_soft_matches","sp_left_sv","sp_right_sv","sp_hom_len","hom"])
    aug = aug.astype(
        {
            "break_pos1": "Int64",
            "split_support": "Int64",
            "soft_support": "Int64",
            "left_soft_matches": "Int64",
            "right_soft_matches": "Int64",
            "sp_left_sv": "Int64",
            "sp_right_sv": "Int64",
            "sp_hom_len": "str",
        }
    )
    out = brk.merge(aug, on="break_pos1", how="left")
    out["b_chr1"] = out["break_chrom1"].apply(lambda x: x[3:]).astype("Int64")
    out = out.sort_values(by=["b_chr1", "break_pos1"])
    out = out.drop(["b_chr1"], axis=1)
    print(f"Augmented SV predictions written to {args.out_table}")
    return out


# -------------------------------------------
# Scaffold pipeline
# -------------------------------------------

aligner = Align.PairwiseAligner(
    mode="local",
    open_gap_score=-10,
    extend_gap_score=-5,
    match_score=2,
    mismatch_score=-10,
)
aligner.wildcard = "?"

def extract_region(fasta, region):
    r = subprocess.run(
        ["samtools", "faidx", fasta, region],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if r.returncode:
        raise RuntimeError(r.stderr)
    return "".join(r.stdout.splitlines()[1:])

# Return list of scaffolds for some breakpoint
def generate_scaffolds(fq1, fq2, out_dir, args):
    shutil.rmtree(out_dir, ignore_errors=True)
    cmd = [
        "python",
        "spades.py",
        "--meta",
        "--pe1-1",
        fq1,
        "--pe1-2",
        fq2,
        "-o",
        out_dir,
        "--pe1-fr",
        "-t",
        str(args.threads)
    ]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode or not os.path.isfile(f"{out_dir}/scaffolds.fasta"):
        raise RuntimeError(f"SPAdes failed:\n{r.stderr, r.stdout}")
    seqs = []
    with open(f"{out_dir}/scaffolds.fasta") as F:
        for block in F.read().split(">")[1:]:
            seqs.append("".join(block.splitlines()[1:]))
    return seqs

# Reformats formatted PairwiseAlignment lists by adding a custom reference chrom and pos, accounting for orientation
def reformat_alignment(arr, ori, ref_chr, ref_start=0, query_start=0):
    output = []
    tgt_fst = None
    tgt_lst = None
    ref_fst = None
    ref_end = None

    for i in range(0, len(arr), 4):
        if i + 2 >= len(arr) or not arr[i].strip():
            continue

        if len(arr[i].split()) < 3:
            continue
        _, tgt_start, seq, *rest_t = arr[i].split()
        tgt_start = int(tgt_start)
        if i == 0:
            query_start = query_start + tgt_start + 1
            tgt_fst = query_start
        seq_len = len(seq) - seq.count("-")
        end_num = query_start + seq_len - 1
        tgt_lst = end_num

        # fixed label + start
        prefix = f"{'target':<15}{query_start} "
        prefix_len = len(prefix)

        # --- parse ref ---
        _, ref_offset, ref_seq, *rest_r = arr[i + 2].split()
        ref_len = len(ref_seq)
        ref_offset = int(ref_offset)
        if i == 0:
            ref_start = (
                (ref_start - 1) + ref_offset - 1000
                if ori
                else (ref_start - 1) + 1000 - ref_offset
            )
            ref_fst = ref_start + 1
        ref_start = ref_start + 1

        ref_end = ref_start + ref_len - ref_seq.count("-") - 1 if ori else ref_start - ref_len + ref_seq.count("-") + 1
        ref_label = f"{'ref':<15}{ref_chr}:{ref_start} "
        ref_row_str = f"{ref_label}{ref_seq} {ref_end}"

        # where the ref bases begin
        ref_seq_col = ref_row_str.index(ref_seq)

        # compute how much to shift the target seq so seq cols match
        pad_spaces = ref_seq_col - prefix_len
        if pad_spaces < 0:
            pad_spaces = 0

        # build target line (label stays flush left)
        first_row_str = prefix + " " * pad_spaces + seq + f" {end_num}"

        # --- build pipes ---
        raw_pipes = arr[i + 1].strip()
        clean_pipes = "".join(ch for ch in raw_pipes if not ch.isdigit())
        # reuse exact same pad that moved the target seq
        pipe_line = " " * (prefix_len + pad_spaces - 1) + clean_pipes

        # --- collect ---
        output.append(first_row_str)
        output.append(pipe_line)
        output.append(ref_row_str)
        output.append("")

        # advance
        query_start = end_num + 1
        ref_start = ref_end if ori else ref_end - 2

    return "\n".join(output), tgt_fst, tgt_lst, ref_fst, ref_end

def test_coordinates(expected_chr, expected_start, expected_end, actual, args):
    region = expected_chr + ":" + str(expected_start) + "-" + str(expected_end)
    expected = extract_region(args.fasta, region).upper()
    if expected != actual:
        print(f"Mismatch\nExpected: {expected}\nActual: {actual}")

def handle_inversion(
    a1_sc_range,
    a1_ref_range,
    a1_is_rev,
    a2_sc_range,
    a2_ref_range,
    a2_is_rev,
    row,
    sc,
    aligner,
):
    seq1 = row["seq1"].upper()
    seq2 = row["seq2"].upper()
    b_pos1 = row["break_pos1"]
    b_pos2 = row["break_pos2"]

    a1_sc_fst, a1_sc_lst = a1_sc_range
    a2_sc_fst, a2_sc_lst = a2_sc_range
    a1_ref_fst, a1_ref_lst = a1_ref_range
    a2_ref_fst, a2_ref_lst = a2_ref_range
    min_sc, is_seq1_min = min(
        (a1_sc_range[0], True), (a2_sc_range[0], False), key=lambda vals: vals[0]
    )
    max_sc, is_seq1_max = max(
        (a1_sc_range[1], True), (a2_sc_range[1], False), key=lambda vals: vals[0]
    )
    min_ref = a1_ref_fst if min_sc == a1_sc_fst else a2_ref_fst
    max_ref = a1_ref_lst if max_sc == a1_sc_lst else a2_ref_lst

    def do_realignment(is_seq1):
        a1, a2 = None, None
        if is_seq1:
            a1 = aligner.align(sc, seq1 if a1_is_rev else rev_comp(seq1))[0]
            a2 = aligner.align(sc, rev_comp(seq2) if a2_is_rev else seq2)[0]
        else:
            a1 = aligner.align(sc, rev_comp(seq1) if a1_is_rev else seq1)[0]
            a2 = aligner.align(sc, seq2 if a2_is_rev else rev_comp(seq2))[0]
        return a1, a2

    # If left side of scaffold is closer to an estimated breakpoint end, then we will realign rev_comp
    if min(abs(min_ref - b_pos1), abs(min_ref - b_pos2)) <= min(
        abs(max_ref - b_pos1), abs(max_ref - b_pos2)
    ):
        a1, a2 = do_realignment(is_seq1_min)
        return a1, a2, is_seq1_min
    else:
        a1, a2 = do_realignment(is_seq1_max)
        return a1, a2, is_seq1_max

def aln_err_density(alns):
    """
    Compute combined density of mismatches ('.') and gaps ('-')
    from pairwise alignments.
    """
    total_errors, total_bases = 0, 0

    for i in range(0, len(alns), 2):
        midline = alns[i + 1].strip()
        total_errors += midline.count(".") + midline.count("-")
        total_bases += len(midline)

    return total_errors / total_bases if total_bases else 0.0

def run_scaffold(args):
    df_all = pd.read_csv(args.file, sep="\t")
    has_homology = ("AA_homology_len" in df_all.columns) and ("AA_homology_seq" in df_all.columns)

    cols = [
        "break_chrom1",
        "break_pos1",
        "break_chrom2",
        "break_pos2",
        "break_sv_type",
        "break_read_support",
        "break_features",
        "break_orientation",
        "amplicon",
    ]
    if has_homology:
        cols.insert(8, "AA_homology_len")
        cols.insert(9, "AA_homology_seq")

    df = df_all[cols].drop_duplicates().reset_index()
    df["ref1"] = (
        df.break_chrom1
        + ":"
        + (df.break_pos1 - 1000).astype(str)
        + "-"
        + (df.break_pos1 + 1000).astype(str)
    )
    df["ref2"] = (
        df.break_chrom2
        + ":"
        + (df.break_pos2 - 1000).astype(str)
        + "-"
        + (df.break_pos2 + 1000).astype(str)
    )
    df["seq1"] = df.ref1.apply(lambda r: extract_region(args.fasta, r))
    df["seq2"] = df.ref2.apply(lambda r: extract_region(args.fasta, r))

    scaffold_log = open(args.scaffold_log + ".txt", "w")
    summary = []

    def fq_names(row):
        b = f"b_{row.break_chrom1}_{row.break_pos1}_{row.break_chrom2}_{row.break_pos2}"
        return f"./fastq/{b}_1.fastq.gz", f"./fastq/{b}_2.fastq.gz"

    for idx, row in df.iterrows():
        fq1, fq2 = fq_names(row)
        bp_out = os.path.join(args.outdir, f"bp_{idx}")
        scaffs = None
        try:
            scaffs = generate_scaffolds(fq1, fq2, bp_out, args)
        except RuntimeError as ex:
            print(
                f"breakpoint {idx + 1}: SpAdes ran into an error. Check that your input arguments are correct. This may also be due to low coverage.\n"
            )
            if args.verbose >= 1:
                print(ex)
            scaffold_log.write(f"=== Breakpoint {idx + 1} ===\n")
            scaffold_log.write(
                "SpAdes ran into an error. Check that your input arguments are correct. This may also be due to low coverage.\n"
            )
            summary.append(
                {
                    "sc_pos1": None,
                    "sc_pos2": None,
                    "sc_hom_len": None,
                    "sc_hom": None,
                }
            )
            continue

        scaffold_log.write(f"=== Breakpoint {idx + 1} ===\n")
        scores = []
        homs = []
        break_ori = row["break_orientation"]
        break_pos1 = row["break_pos1"]
        break_pos2 = row["break_pos2"]
        for i, sc in enumerate(scaffs):
            a1_pos = aligner.align(sc, row.seq1.upper())
            a1_neg = aligner.align(sc, rev_comp(row.seq1.upper()))
            a2_pos = aligner.align(sc, row.seq2.upper())
            a2_neg = aligner.align(sc, rev_comp(row.seq2.upper()))

            if (len(a1_pos) == 0 or len(a1_neg) == 0) or (len(a2_pos) == 0 or len(a2_neg) == 0):
                print("Alignment score too low")
                scaffold_log.write(
                f"-- scaffold {i + 1} --\nscaffold {sc}\nAlignment score too low\n"
                )
                continue

            a1_pos = a1_pos[0]
            a1_neg = a1_neg[0]
            a2_pos = a2_pos[0]
            a2_neg = a2_neg[0]

            a1, a1_is_rev = max(
                (a1_pos, True), (a1_neg, False), key=lambda aln: aln[0].score
            )
            a2, a2_is_rev = max(
                (a2_pos, True), (a2_neg, False), key=lambda aln: aln[0].score
            )
            a1_out, a1_fst, a1_lst, a1_ref_fst, a1_ref_lst = reformat_alignment(
                format(a1).split("\n"),
                (a1_pos.score >= a1_neg.score),
                row["break_chrom1"],
                row["break_pos1"],
            )
            a2_out, a2_fst, a2_lst, a2_ref_fst, a2_ref_lst = reformat_alignment(
                format(a2).split("\n"),
                (a2_pos.score >= a2_neg.score),
                row["break_chrom2"],
                row["break_pos2"],
            )

            hom = ""
            hom_len = None
            l_ref_pos = None
            r_ref_pos = None
            is_pos_valid = True
            if (
                (row["break_orientation"] == "++" or row["break_orientation"] == "--")
                and row["break_sv_type"] != "interchromosomal"
                and (
                    (a1_fst >= a2_fst and a1_lst <= a2_lst)
                    or (a2_fst >= a1_fst and a2_lst <= a1_lst)
                )
            ):
                print("Inversion detected.\n")
                a1, a2, is_seq1 = handle_inversion(
                    (a1_fst, a1_lst),
                    (a1_ref_fst, a1_ref_lst),
                    a1_is_rev,
                    (a2_fst, a2_lst),
                    (a2_ref_fst, a2_ref_lst),
                    a2_is_rev,
                    row,
                    sc,
                    aligner,
                )
                a1_out, a1_fst, a1_lst, a1_ref_fst, a1_ref_lst = reformat_alignment(
                    format(a1).split("\n"),
                    (is_seq1 and a1_is_rev) or (not is_seq1 and not a1_is_rev),
                    row["break_chrom1"],
                    row["break_pos1"],
                )
                a2_out, a2_fst, a2_lst, a2_ref_fst, a2_ref_lst = reformat_alignment(
                    format(a2).split("\n"),
                    (is_seq1 and not a2_is_rev) or (not is_seq1 and a2_is_rev),
                    row["break_chrom2"],
                    row["break_pos2"],
                )
                print("Inversion handled.\n")

            if (
                a1.score <= 50
                or a2.score <= 50
                or aln_err_density(a1_out.split("\n")) >= 0.1
                or aln_err_density(a2_out.split("\n")) >= 0.1
            ):
                hom = None
            elif (a1_fst >= a2_fst and a1_lst <= a2_lst) or (
                a2_fst >= a1_fst and a2_lst <= a1_lst
            ):
                hom = None
                print(
                    f"Inversion may not have been correctly detected for scaffold {i + 1}.\n"
                )
            elif a1_lst + 1 == a2_fst:
                hom = "N/A"
                hom_len = 0
                l_ref_pos = a1_ref_lst
                r_ref_pos = a2_ref_fst
            elif a2_lst + 1 == a1_fst:
                hom = "N/A"
                hom_len = 0
                l_ref_pos = a1_ref_fst
                r_ref_pos = a2_ref_lst
            elif a1_lst < a2_fst:
                hom = sc[a1_lst : a2_fst - 1]
                hom_len = -len(hom)
                l_ref_pos = a1_ref_lst
                r_ref_pos = a2_ref_fst
            elif a2_lst < a1_fst:
                hom = sc[a2_lst : a1_fst - 1]
                hom_len = -len(hom)
                l_ref_pos = a1_ref_fst
                r_ref_pos = a2_ref_lst
            elif a1_lst <= a2_lst and a1_lst >= a2_fst:
                hom = sc[a2_fst - 1 : a1_lst]
                hom_len = len(hom)
                l_ref_pos = a1_ref_lst
                r_ref_pos = a2_ref_fst
            else:
                hom = sc[a1_fst - 1 : a2_lst]
                hom_len = len(hom)
                l_ref_pos = a1_ref_fst
                r_ref_pos = a2_ref_lst
            
            if hom is not None:
                left_within_20 = (l_ref_pos >= break_pos1-20) & (l_ref_pos <= break_pos1+20)
                right_within_20 = (r_ref_pos >= break_pos2-20) & (r_ref_pos <= break_pos2+20)
                left_match_ori = (l_ref_pos <= break_pos1) if break_ori[0] == "-" else (l_ref_pos >= break_pos1)
                right_match_ori = (r_ref_pos <= break_pos2) if break_ori[1] == "-" else (r_ref_pos >= break_pos2)
                is_pos_valid = (left_within_20 or left_match_ori) and (right_within_20 or right_match_ori)

            scaffold_log.write(
                f"-- scaffold {i + 1} --\nscaffold {sc}\nscaffold length {len(sc)} bp\n\nref1->scaffold:\n{a1_out}\n"
            )
            scaffold_log.write(f"ref2->scaffold:\n{a2_out}\n")
            scaffold_log.write(f"Homology/Insertion: {hom}\n")
            scaffold_log.write(f"Homology Length: {hom_len}\n")
            scores.append((a1.score, a2.score))

            if hom is not None and hom != "N/A" and a1_lst >= a2_lst:
                hom = rev_comp(hom)

            if is_pos_valid:
                homs.append((hom_len, hom, l_ref_pos, r_ref_pos))
            else:
                homs.append((None, None, None, None))

        if len(scores) != 0:
            li = max(range(len(scores)), key=lambda j: scores[j][0])
            ri = max(range(len(scores)), key=lambda j: scores[j][1])
            print(
                f"breakpoint {idx + 1}: left={li}({scores[li]}), right={ri}({scores[ri]})"
            )

        best_prediction = next((s for s in homs if s[1]), ("", "", "", ""))
        summary.append(
            {
                "sc_pos1": best_prediction[2],
                "sc_pos2": best_prediction[3],
                "sc_hom_len": best_prediction[0],
                "sc_hom": best_prediction[1],
            }
        )
        homs = []

    scaffold_log.close()

    base_cols = [
        "break_chrom1",
        "break_pos1",
        "break_chrom2",
        "break_pos2",
        "break_sv_type",
        "break_read_support",
        "break_features",
        "break_orientation",
        "amplicon",
    ]
    if has_homology:
        base_cols.insert(8, "AA_homology_len")
        base_cols.insert(9, "AA_homology_seq")

    base = df[base_cols]
    aug = pd.DataFrame(summary)
    out = pd.concat([base.reset_index(drop=True), aug], axis=1)
    out["sc_hom_len"] = out["sc_hom_len"].replace("", np.nan)
    out["sc_hom_len"] = out["sc_hom_len"].astype(float).astype("Int64")
    out["b_chr1"] = out["break_chrom1"].apply(lambda x: x[3:]).astype("Int64")
    out = out.sort_values(by=["b_chr1", "break_pos1"])
    out = out.drop(["b_chr1"], axis=1)
    print(f"Augmented scaffold predictions written to {args.out_table}")
    return out

def main():
    p = argparse.ArgumentParser(
        description="SV split read analysis OR scaffold reconstruction"
    )
    p.add_argument("file", help="input TSV")
    p.add_argument(
        "--mode",
        choices=["split", "scaffold", "both"],
        default="split",
        help="which pipeline to run",
    )
    p.add_argument(
        "--out-table",
        default="augmented_predictions.tsv",
        help="path for augmented output TSV",
    )
    p.add_argument(
        "--split-log", default="split_read_alignments", help="refine pipeline log"
    )
    p.add_argument(
        "--scaffold-log",
        default="scaffold_alignments",
        help="scaffold pipeline log",
    )
    p.add_argument("--outdir", default="out", help="SPAdes output base directory")
    p.add_argument("-l", "--list", action="store_true", help="list breakpoints")
    p.add_argument(
        "-b",
        "--breakpoints",
        nargs="+",
        type=int,
        help="which breakpoint indices to process",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="verbosity (use -vv for more)",
    )
    p.add_argument(
        "-t",
        "--threads",
        type=int,
        default=16,
        help="how many threads to use for the assembler",
    )
    p.add_argument("--fasta", help="indexed FASTA for extract_region")

    args = p.parse_args()
    args.out_table = args.out_table + ".tsv"
    if args.mode == "split":
        out = run_split(args)
        if out is not None:
            out.to_csv(args.out_table, sep="\t", index=False)
    elif args.mode == "scaffold":
        if not args.fasta:
            p.error("--fasta is required in scaffold mode")
        out = run_scaffold(args)
        out.to_csv(args.out_table, sep="\t", index=False)
    else:
        if not args.fasta:
            p.error("--fasta is required in scaffold mode")
        out_sc = run_scaffold(args)
        out_split = run_split(args)
        final = pd.concat(
            [
                out_sc.reset_index(drop=True),
                out_split[
                    [
                        "split_support",
                        "soft_support",
                        "left_soft_matches",
                        "right_soft_matches",
                        "sp_left_sv",
                        "sp_right_sv",
                        "sp_hom_len",
                        "hom",
                    ]
                ].reset_index(drop=True),
            ],
            axis=1,
        )
        final.to_csv(args.out_table, sep="\t", index=False)

if __name__ == "__main__":
    main()
