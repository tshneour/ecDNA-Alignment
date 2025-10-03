#!/usr/bin/env python3
import os
import argparse
import warnings
import pandas as pd
from Bio import Align
import subprocess, shutil
from natsort import index_natsorted
import numpy as np

# suppress pandas warnings
warnings.filterwarnings("ignore")
pd.set_option("mode.chained_assignment", None)

# -------------------------------------------
# Utility functions
# -------------------------------------------


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
    "query_short",
    "query_chrom",
    "query_pos",
    "query_cigar",
    "hom_clip_match",
    "ins_clip_match",
    "rev_clipped",
    "query_aln_sub",
]


def refine_step1(reads, all_reads, verbose):
    reads = reads[reads["query_cigar"].str.contains(r"[SH]")].copy()
    reads["begin"] = reads["query_cigar"].str.contains(r"^\d+[SH]", regex=True)
    reads["end"] = reads["query_cigar"].str.contains(r"\d+[SH]$", regex=True)
    reads["sv_end"] = reads["query_pos"].where(reads["begin"])
    reads["sv_end"].fillna(reads["query_end"].where(reads["end"]), inplace=True)
    reads.dropna(subset=["sv_end"], inplace=True)

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

    return best, reads


def check_overlap(left, right):
    # filter splits
    left = left[
        ~left["split"]
        | (
            left["split"]
            & left["query_name"].isin(right.loc[right["split"], "query_name"])
        )
    ]
    right = right[
        ~right["split"]
        | (
            right["split"]
            & right["query_name"].isin(left.loc[left["split"], "query_name"])
        )
    ]

    # compute clipped sequences
    for df in (left, right):
        df["clip_len"] = df["query_aln_full"].str.len() - df["query_aln_sub"].str.len()
        # print(df.dtypes)
        # print(type(df["query_aln_full"]))
        # print(df["query_aln_full"])
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

    results = []
    for lv, ldf in left.groupby("sv_end"):
        for rv, rdf in right.groupby("sv_end"):
            # pick the longest sub-alignment from each
            fr = ldf.loc[ldf["query_aln_sub"].str.len().idxmax()]
            lr = rdf.loc[rdf["query_aln_sub"].str.len().idxmax()]

            # CIGAR strings already in pos dir so just use soft clipped side to find sv ori
            l_ori = "+" if fr["end"] else "-"
            r_ori = "+" if lr["end"] else "-"
            sv_ori = l_ori + r_ori
            if sv_ori != fr["break_orientation"]:
                continue

            seq1 = fr["query_aln_sub"] if fr["end"] else rev_comp(fr["query_aln_sub"])
            seq2 = lr["query_aln_sub"] if lr["begin"] else rev_comp(lr["query_aln_sub"])

            hom = get_homology(seq1, seq2)
            # if hom == "":
            #     hom = "N/A"
            hlen = len(hom)
            s1_no = seq1[:-hlen]
            s2_no = seq2[hlen:]

            # homology matching
            for df, target_no, name in (
                (ldf, s2_no, "hom_clip_match"),
                (rdf, s1_no, "hom_clip_match"),
            ):
                df["rev_clipped"] = df["clipped"].where(
                    df["end"] if df is ldf else df["begin"], rev_comp
                )
                df[name] = df["rev_clipped"].apply(lambda x: target_no.startswith(x)) if df is ldf else df["rev_clipped"].apply(lambda x: target_no.endswith(x))
                # fix split reads if necessary
                df.loc[df["split"], name] = df.loc[df["split"]].apply(
                    lambda x: (
                        x["query_short"]
                        in (
                            (rdf if df is ldf else ldf).loc[
                                (rdf if df is ldf else ldf)["split"], "query_short"
                            ]
                        ).values
                        if "H" in x["query_cigar"]
                        else x[name]
                    ),
                    axis=1,
                )

            hsum_l = ldf["hom_clip_match"].sum()
            hsum_r = rdf["hom_clip_match"].sum()

            # left_clip     CTGTCACTGACTGAATTAATTTTACCTCATATGTCTGG
            # right_clip    CTGTCACTGACTGAATTAATTTTACCTCATATGTC
            # seq1          GTGCCTATGTGTTTATCTCCACAGTGCAATGTATTGGTTACATAATAGCTGTCACTGACTGAATTAATTTTACCTCATATGTCTGG
            # seq2          TGGGCAAAATCTTAAGCTCAGTTTTTTAACTTTATTTTTGGTCCATCTTGATTATAATTCTATTT
            # ins           CTGTCACTGACTGAATTAATTTTACCTCATATGTC

            # insertion matching
            c1 = ldf.loc[ldf["rev_clipped"].str.len().idxmax(), "rev_clipped"]
            c2 = rdf.loc[rdf["rev_clipped"].str.len().idxmax(), "rev_clipped"]
            ins = get_homology(c2, c1)
            ilen = len(ins)
            # print(c1)
            # print(c2)
            # print(seq1)
            # print(seq2)
            # print(ins, "\n")
            # print(ldf["rev_clipped"].str.slice(start=ilen))
            # print(rdf["rev_clipped"].str.slice(stop=-ilen))

            ldf["ins_clip_match"] = (
                ldf["rev_clipped"]
                .str.slice(start=ilen)
                .apply(lambda x: seq2.startswith(x) if len(x) != 0 else True)
            )
            rdf["ins_clip_match"] = (
                rdf["rev_clipped"]
                .str.slice(stop=-ilen)
                .apply(lambda x: (seq1.endswith(x) if len(x) != 0 else True) if ilen != 0 else False)
            )
            # print(ldf["ins_clip_match"])
            # print(rdf["ins_clip_match"])
            # fix split reads for insertion
            for df in (ldf, rdf):
                df.loc[df["split"], "ins_clip_match"] = df.loc[df["split"]].apply(
                    lambda x: (
                        x["query_short"]
                        in (
                            (rdf if df is ldf else ldf).loc[
                                (rdf if df is ldf else ldf)["split"], "query_short"
                            ]
                        ).values
                        if "H" in x["query_cigar"]
                        else x["ins_clip_match"]
                    ),
                    axis=1,
                )

            isum_l = ldf["ins_clip_match"].sum()
            isum_r = rdf["ins_clip_match"].sum()

            # Output collection
            total_reads = len(ldf) + len(rdf)
            num_split_reads = (ldf["split"].sum() + rdf["split"].sum()) / 2

            left_disc_mask = (ldf["proper_pair"] == "Discordant") & (~ldf["split"])
            right_disc_mask = (rdf["proper_pair"] == "Discordant") & (~rdf["split"])
            left_conc_mask = (ldf["proper_pair"] == "Concordant") & (~ldf["split"])
            right_conc_mask = (rdf["proper_pair"] == "Concordant") & (~rdf["split"])
            left_disc_matches_hom = ldf.drop_duplicates().loc[left_disc_mask, "hom_clip_match"].sum()
            right_disc_matches_hom = rdf.drop_duplicates().loc[right_disc_mask, "hom_clip_match"].sum()
            left_disc_matches_ins = ldf.drop_duplicates().loc[left_disc_mask, "ins_clip_match"].sum()
            right_disc_matches_ins = rdf.drop_duplicates().loc[right_disc_mask, "ins_clip_match"].sum()
            left_conc_matches_hom = ldf.loc[left_conc_mask, "hom_clip_match"].sum()
            right_conc_matches_hom = rdf.loc[right_conc_mask, "hom_clip_match"].sum()
            left_conc_matches_ins = ldf.loc[left_conc_mask, "ins_clip_match"].sum()
            right_conc_matches_ins = rdf.loc[right_conc_mask, "ins_clip_match"].sum()
            
            num_left_soft_reads = len(ldf) - ldf["split"].sum() - (len(ldf.drop_duplicates().loc[left_disc_mask]))
            num_right_soft_reads = len(rdf) - rdf["split"].sum() - (len(rdf.drop_duplicates().loc[right_disc_mask]))

            # Filter out cases with poor soft clip matches
            if num_split_reads == 0:
                hom_left = hsum_l / num_left_soft_reads
                hom_right = hsum_r / num_right_soft_reads
                ins_left = isum_l / num_left_soft_reads
                ins_right = isum_r / num_right_soft_reads
                if (hom_left < 0.33 and hom_right < 0.33) and (
                    ins_left < 0.33 and ins_right < 0.33
                ):
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
                        "left_len": [len(ldf)],
                        "right_len": [len(rdf)],
                        "split_support": [num_split_reads // 2],
                        "soft_support": [num_left_soft_reads + num_right_soft_reads],
                        "left_soft_matches_hom": [left_disc_matches_hom + left_conc_matches_hom],
                        "right_soft_matches_hom": [right_disc_matches_hom + right_conc_matches_hom],
                        "left_soft_matches_ins": [left_disc_matches_ins + left_conc_matches_ins],
                        "right_soft_matches_ins": [right_disc_matches_ins + right_conc_matches_ins],
                        "total_reads": [total_reads],
                        # <-- store DataFrame directly, not as tuple
                        "left": [ldf.copy()],
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
    all_reads = pd.read_csv(args.file, sep="\t")
    svs = all_reads.groupby("break_pos1")

    split_log = open(args.split_log + ".txt", "w")
    summary = []

    if args.list:
        tbl = (
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

        # Inversion case
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

        _, lgrp = refine_step1(left, all_reads, args.verbose)
        _, rgrp = refine_step1(right, all_reads, args.verbose)

        res = check_overlap(lgrp, rgrp)

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
        res = res.sort_values(["split_support", "total_%", "total_reads"], ascending=False)

        res = res.head(3)

        if res is None or len(res) == 0:
            print("No overlap found")
            continue

        # res["homology"] = res["homology"].apply(lambda x: x if break_ori[0] == "+" else rev_comp(x))
        # res["insertion"] = res["insertion"].apply(lambda x: x if break_ori[0] == "+" else rev_comp(x))

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

        for i in range(len(res) - 1, -1, -1):
            row = res.iloc[i]
            print("_" * 80)
            print(row.drop(["left", "right"]).to_string())
            print("Left reads:\n", row["left"][print_cols2].to_string())
            print("Right reads:\n", row["right"][print_cols2].to_string())
            print("_" * 80, "\n")

        top = res.iloc[0]
        # print(top.index)
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
                else -len(top.insertion),
                "hom": top.homology if top["hom_%"] > top["ins_%"] else top.insertion,
            }
        )

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

    brk = (
        svs[
            [
                "break_chrom1",
                "break_pos1",
                "break_chrom2",
                "break_pos2",
                "break_sv_type",
                "break_read_support",
                "break_features",
                "break_orientation",
                "AA_homology_len",
                "AA_homology_seq",
                "amplicon",
            ]
        ]
        .first()
        .reset_index(drop=True)
    )
    aug = pd.DataFrame(summary)
    aug = aug.astype(
        {
            "break_pos1": "Int64",
            "split_support": "Int64",
            "soft_support": "Int64",
            "left_soft_matches": "Int64",
            "right_soft_matches": "Int64",
            "sp_left_sv": "Int64",
            "sp_right_sv": "Int64",
            "sp_hom_len": "Int64",
        }
    )
    # aug = aug[["break_chrom1", "break_pos1", "break_chrom2", "break_pos2", "break_sv_type", "break_read_support", "break_features", "break_orientation"]]
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
def generate_scaffolds(fq1, fq2, out_dir):
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
        "--pe1-fr"
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

        # --- parse target ---
        # if i == 0:
        #     print('\n'.join(arr))
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
                (ref_start - 1) + ref_offset - 500
                if ori
                else (ref_start - 1) + 500 - ref_offset
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
            # Not 100% about this
            a1 = aligner.align(sc, rev_comp(seq1) if a1_is_rev else seq1)[0]
            a2 = aligner.align(sc, seq2 if a2_is_rev else rev_comp(seq2))[0]
        return a1, a2

    # If left side of scaffold is closer to an estimated breakpoint end, then we will realign rev_comp
    if min(abs(min_ref - b_pos1), abs(min_ref - b_pos2)) <= min(
        abs(max_ref - b_pos1), abs(max_ref - b_pos2)
    ):
        a1, a2 = do_realignment(is_seq1_min)
        # print(f"{format(a1).split("\n")}\n{format(a2).split("\n")}\n")
        # print(f"{format(a1).split("\n")}\n{format(a2).split("\n")}\n")
        return a1, a2, is_seq1_min
    else:
        a1, a2 = do_realignment(is_seq1_max)
        # print(f"{format(a1).split("\n")}\n{format(a2).split("\n")}\n")
        # print(f"{format(a1).split("\n")}\n{format(a2).split("\n")}\n")
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
    df = (
        pd.read_csv(args.file, sep="\t")[
            [
                "break_chrom1",
                "break_pos1",
                "break_chrom2",
                "break_pos2",
                "break_sv_type",
                "break_read_support",
                "break_features",
                "break_orientation",
                "AA_homology_len",
                "AA_homology_seq",
                "amplicon",
            ]
        ]
        .drop_duplicates()
        .reset_index()
    )
    df["ref1"] = (
        df.break_chrom1
        + ":"
        + (df.break_pos1 - 500).astype(str)
        + "-"
        + (df.break_pos1 + 500).astype(str)
    )
    df["ref2"] = (
        df.break_chrom2
        + ":"
        + (df.break_pos2 - 500).astype(str)
        + "-"
        + (df.break_pos2 + 500).astype(str)
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
            scaffs = generate_scaffolds(fq1, fq2, bp_out)
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
            # print(sc)
            a1_pos = aligner.align(sc, row.seq1.upper())[0]
            a1_neg = aligner.align(sc, rev_comp(row.seq1.upper()))[0]
            a2_pos = aligner.align(sc, row.seq2.upper())[0]
            a2_neg = aligner.align(sc, rev_comp(row.seq2.upper()))[0]
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
            first_aln_len = len(a1_out.split("\n")[2].split()[2]) - 1
            # test_coordinates(
            #     df["break_chrom1"][idx],
            #     a1_ref_fst - first_aln_len
            #     if (a1_pos.score < a1_neg.score)
            #     else a1_ref_fst,
            #     a1_ref_fst
            #     if (a1_pos.score < a1_neg.score)
            #     else a1_ref_fst + first_aln_len,
            #     rev_comp(a1_out.split("\n")[2].split()[2])
            #     if (a1_pos.score < a1_neg.score)
            #     else a1_out.split("\n")[2].split()[2],
            #     args,
            # )
            a2_out, a2_fst, a2_lst, a2_ref_fst, a2_ref_lst = reformat_alignment(
                format(a2).split("\n"),
                (a2_pos.score >= a2_neg.score),
                row["break_chrom2"],
                row["break_pos2"],
            )
            second_aln_len = len(a2_out.split("\n")[2].split()[2]) - 1
            # test_coordinates(
            #     df["break_chrom2"][idx],
            #     a2_ref_fst - second_aln_len
            #     if (a2_pos.score < a2_neg.score)
            #     else a2_ref_fst,
            #     a2_ref_fst
            #     if (a2_pos.score < a2_neg.score)
            #     else a2_ref_fst + second_aln_len,
            #     rev_comp(a2_out.split("\n")[2].split()[2])
            #     if (a2_pos.score < a2_neg.score)
            #     else a2_out.split("\n")[2].split()[2],
            #     args,
            # )

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
                first_aln_len = len(a1_out.split("\n")[2].split()[2]) - 1
                # test_coordinates(
                #     df["break_chrom1"][idx],
                #     a1_ref_fst - first_aln_len
                #     if not ((is_seq1 and a1_is_rev) or (not is_seq1 and not a1_is_rev))
                #     else a1_ref_fst,
                #     a1_ref_fst
                #     if not ((is_seq1 and a1_is_rev) or (not is_seq1 and not a1_is_rev))
                #     else a1_ref_fst + first_aln_len,
                #     rev_comp(a1_out.split("\n")[2].split()[2])
                #     if not ((is_seq1 and a1_is_rev) or (not is_seq1 and not a1_is_rev))
                #     else a1_out.split("\n")[2].split()[2],
                #     args,
                # )
                a2_out, a2_fst, a2_lst, a2_ref_fst, a2_ref_lst = reformat_alignment(
                    format(a2).split("\n"),
                    (is_seq1 and not a2_is_rev) or (not is_seq1 and a2_is_rev),
                    row["break_chrom2"],
                    row["break_pos2"],
                )
                second_aln_len = len(a2_out.split("\n")[2].split()[2]) - 1
                # test_coordinates(
                #     df["break_chrom2"][idx],
                #     a2_ref_fst - second_aln_len
                #     if not ((is_seq1 and not a2_is_rev) or (not is_seq1 and a2_is_rev))
                #     else a2_ref_fst,
                #     a2_ref_fst
                #     if not ((is_seq1 and not a2_is_rev) or (not is_seq1 and a2_is_rev))
                #     else a2_ref_fst + second_aln_len,
                #     rev_comp(a2_out.split("\n")[2].split()[2])
                #     if not ((is_seq1 and not a2_is_rev) or (not is_seq1 and a2_is_rev))
                #     else a2_out.split("\n")[2].split()[2],
                #     args,
                # )
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

        li = max(range(len(scores)), key=lambda j: scores[j][0])
        ri = max(range(len(scores)), key=lambda j: scores[j][1])
        # summary.append(
        #     {
        #         "best_left_scaffold": li,
        #         "best_left_scores": scores[li],
        #         "best_right_scaffold": ri,
        #         "best_right_scores": scores[ri],
        #     }
        # )
        # print([s for s in homs if s[1]])
        best_prediction = next((s for s in homs if s[1]), ("", "", "", ""))
        # print(best_prediction)
        summary.append(
            {
                "sc_pos1": best_prediction[2],
                "sc_pos2": best_prediction[3],
                "sc_hom_len": best_prediction[0],
                "sc_hom": best_prediction[1],
            }
        )
        print(
            f"breakpoint {idx + 1}: left={li}({scores[li]}), right={ri}({scores[ri]})"
        )
        homs = []

    scaffold_log.close()

    base = df.drop(columns=["ref1", "ref2", "seq1", "seq2"])
    base = df[
        [
            "break_chrom1",
            "break_pos1",
            "break_chrom2",
            "break_pos2",
            "break_sv_type",
            "break_read_support",
            "break_features",
            "break_orientation",
            "AA_homology_len",
            "AA_homology_seq",
            "amplicon",
        ]
    ]
    aug = pd.DataFrame(summary)
    out = pd.concat([base.reset_index(drop=True), aug], axis=1)
    out["sc_hom_len"] = out["sc_hom_len"].replace("", np.nan)
    out["sc_hom_len"] = out["sc_hom_len"].astype(float).astype("Int64")
    out["b_chr1"] = out["break_chrom1"].apply(lambda x: x[3:]).astype("Int64")
    out = out.sort_values(by=["b_chr1", "break_pos1"])
    out = out.drop(["b_chr1"], axis=1)
    print(f"Augmented scaffold predictions written to {args.out_table}")
    return out


# -------------------------------------------
# Main & CLI
# -------------------------------------------


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
    # refine flags
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
    # scaffold flags
    p.add_argument("--fasta", help="indexed FASTA for extract_region")

    args = p.parse_args()
    args.out_table = args.out_table + ".tsv"
    if args.mode == "split":
        out = run_split(args)
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
