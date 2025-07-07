#!/usr/bin/env python3
import os
import argparse
import warnings
import pandas as pd
from Bio import Align
import subprocess, shutil

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
    "query_name","proper_pair","query_pos","query_end",
    "query_cigar","split","query_aln_full",
]
print_cols2 = [
    "query_short","query_chrom","query_pos","query_cigar",
    "hom_clip_match","ins_clip_match","rev_clipped","query_aln_sub",
]

def refine_step1(reads, all_reads, verbose):
    reads = reads[reads["query_cigar"].str.contains(r"[SH]")].copy()
    reads["begin"] = reads["query_cigar"].str.contains(r"^\d+[SH]", regex=True)
    reads["end"]   = reads["query_cigar"].str.contains(r"\d+[SH]$", regex=True)
    reads["sv_end"] = reads["query_pos"].where(reads["begin"])
    reads["sv_end"].fillna(reads["query_end"].where(reads["end"]), inplace=True)
    reads.dropna(subset=["sv_end"], inplace=True)

    if verbose >= 2:
        print(reads.sort_values("sv_end")[["sv_end",*print_cols]].to_string(index=False),"\n")

    cands = reads["sv_end"].value_counts()
    best = cands.index[0] if len(cands)>0 else None

    if verbose == 1 and best is not None:
        sel = reads[reads["sv_end"]==best][["sv_end",*print_cols]]
        disc = sel[sel["proper_pair"]=="Discordant"]
        mates = all_reads[
            all_reads["query_name"].isin(disc.query_name)
            & ~all_reads["query_aln_full"].isin(disc.query_aln_full)
        ][print_cols]
        print(sel.to_string(index=False),"\nMates:\n",mates.to_string(index=False),"\n")

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
        df["clipped"] = df.apply(
            lambda x: (
                x["query_aln_full"][: x["clip_len"]]
                if x["begin"]
                else x["query_aln_full"][-x["clip_len"] or len(x["query_aln_full"]) :]
            ),
            axis=1,
        )

    results = []
    for lv, ldf in left.groupby("sv_end"):
        for rv, rdf in right.groupby("sv_end"):
            # pick the longest sub-alignment from each
            fr = ldf.loc[ldf["query_aln_sub"].str.len().idxmax()]
            lr = rdf.loc[rdf["query_aln_sub"].str.len().idxmax()]

            seq1 = fr["query_aln_sub"] if fr["end"] else rev_comp(fr["query_aln_sub"])
            seq2 = lr["query_aln_sub"] if lr["begin"] else rev_comp(lr["query_aln_sub"])

            hom = get_homology(seq1, seq2)
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
                df[name] = df["rev_clipped"].apply(lambda x: target_no.startswith(x))
                # fix split reads if necessary
                df.loc[df["split"], name] = df.loc[df["split"]].apply(
                    lambda x: (
                        x["query_short"]
                        in (
                            (rdf if df is ldf else ldf)
                            .loc[(rdf if df is ldf else ldf)["split"], "query_short"]
                        ).values
                        if "H" in x["query_cigar"]
                        else x[name]
                    ),
                    axis=1,
                )

            hsum_l = ldf["hom_clip_match"].sum()
            hsum_r = rdf["hom_clip_match"].sum()

            # insertion matching
            c1 = ldf.loc[ldf["rev_clipped"].str.len().idxmax(), "rev_clipped"]
            c2 = rdf.loc[rdf["rev_clipped"].str.len().idxmax(), "rev_clipped"]
            ins = get_homology(c2, c1)
            ilen = len(ins)

            ldf["ins_clip_match"] = ldf["rev_clipped"].str.slice(start=ilen).apply(
                lambda x: seq2.startswith(x)
            )
            rdf["ins_clip_match"] = rdf["rev_clipped"].str.slice(stop=-ilen).apply(
                lambda x: seq1.endswith(x)
            )
            # fix split reads for insertion
            for df in (ldf, rdf):
                df.loc[df["split"], "ins_clip_match"] = df.loc[df["split"]].apply(
                    lambda x: (
                        x["query_short"]
                        in (
                            (rdf if df is ldf else ldf)
                            .loc[(rdf if df is ldf else ldf)["split"], "query_short"]
                        ).values
                        if "H" in x["query_cigar"]
                        else x["ins_clip_match"]
                    ),
                    axis=1,
                )

            isum_l = ldf["ins_clip_match"].sum()
            isum_r = rdf["ins_clip_match"].sum()

            total_reads = len(ldf) + len(rdf)
            split_matches = ldf["split"].sum() + rdf["split"].sum()

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
                        "split_matches": [split_matches],
                        "total_reads": [total_reads],
                        # <-- store DataFrame directly, not as tuple
                        "left": [ldf.copy()],
                        "right": [rdf.copy()],
                    }
                )
            )

    out = pd.concat(results, ignore_index=True)
    out["hom_%"] = (out["hom_sum_left"] + out["hom_sum_right"]) / (
        out["left_len"] + out["right_len"]
    )
    out["ins_%"] = (out["ins_sum_left"] + out["ins_sum_right"]) / (
        out["left_len"] + out["right_len"]
    )
    out["total_%"] = out[["hom_%", "ins_%"]].max(axis=1)
    return out.sort_values(["total_%", "split_matches", "total_reads"], ascending=False)

def run_refine(args):
    all_reads = pd.read_csv(args.file, sep="\t")
    svs = all_reads.groupby("break_pos1")

    split_log = open(args.split_log, "w")
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

        left = sv[
            (sv["query_pos"] >= bp - 350) & (sv["query_end"] <= bp + 350)
        ]
        right = sv[
            (sv["query_pos"] >= sv["break_pos2"].iloc[0] - 350)
            & (sv["query_end"] <= sv["break_pos2"].iloc[0] + 350)
        ]

        _, lgrp = refine_step1(left, all_reads, args.verbose)
        _, rgrp = refine_step1(right, all_reads, args.verbose)

        res = check_overlap(lgrp, rgrp).head(3)

        split_log.write(f"=== Breakpoint {bp} ===\n")
        for idx, row in res.iterrows():
            split_log.write(f"-- Candidate {idx} --\n")
            split_log.write(row.drop(["left", "right"]).to_string() + "\n")
            split_log.write("Left reads:\n" + row["left"][print_cols2].to_string() + "\n")
            split_log.write("Right reads:\n" + row["right"][print_cols2].to_string() + "\n\n")

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
                "left_sv": int(top.left_sv),
                "right_sv": int(top.right_sv),
                "homology": top.homology,
                "insertion": top.insertion,
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
                "AA_homology_seq",
            ]
        ]
        .first()
        .reset_index()
    )
    aug = pd.DataFrame(summary)
    out = brk.merge(aug, on="break_pos1", how="left")
    out.to_csv(args.out_table, sep="\t", index=False)
    print(f"Augmented SV predictions written to {args.out_table}")

# -------------------------------------------
# Scaffold pipeline
# -------------------------------------------

aligner = Align.PairwiseAligner(
    mode="local",
    open_gap_score=-10,
    extend_gap_score=-5,
    match_score=2,
    mismatch_score=-2,
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

def generate_scaffolds(fq1, fq2, out_dir):
    shutil.rmtree(out_dir, ignore_errors=True)
    cmd = [
        "python",
        "./SPAdes-4.2.0-Linux/bin/spades.py",
        "--meta",
        "--pe1-1",
        fq1,
        "--pe1-2",
        fq2,
        "-o",
        out_dir,
    ]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode or not os.path.isfile(f"{out_dir}/scaffolds.fasta"):
        raise RuntimeError(f"SPAdes failed:\n{r.stderr, r.stdout}")
    seqs = []
    with open(f"{out_dir}/scaffolds.fasta") as F:
        for block in F.read().split(">")[1:]:
            seqs.append("".join(block.splitlines()[1:]))
    return seqs

def run_scaffold(args):
    df = pd.read_csv(args.file, sep="\t")[["break_chrom1", "break_pos1", "break_chrom2", "break_pos2"]].drop_duplicates().reset_index()
    df["ref1"] = df.break_chrom1 + ":" + (df.break_pos1 - 350).astype(str) + "-" + (
        df.break_pos1 + 350
    ).astype(str)
    df["ref2"] = df.break_chrom2 + ":" + (df.break_pos2 - 350).astype(str) + "-" + (
        df.break_pos2 + 350
    ).astype(str)
    df["seq1"] = df.ref1.apply(lambda r: extract_region(args.fasta, r))
    df["seq2"] = df.ref2.apply(lambda r: extract_region(args.fasta, r))

    scaffold_log = open(args.scaffold_log, "w")
    summary = []

    def fq_names(row):
        b = f"b_{row.break_chrom1}:{row.break_pos1}_{row.break_chrom2}:{row.break_pos2}"
        return f"./fastq/{b}_1.fastq.gz", f"./fastq/{b}_2.fastq.gz"

    for idx, row in df.iterrows():
        fq1, fq2 = fq_names(row)
        bp_out = os.path.join(args.outdir, f"bp_{idx}")
        scaffs = generate_scaffolds(fq1, fq2, bp_out)

        scaffold_log.write(f"=== Breakpoint {idx} ===\n")
        scores = []
        for i, sc in enumerate(scaffs):
            a1 = aligner.align(row.seq1.upper(), sc)[0]
            a2 = aligner.align(row.seq2.upper(), sc)[0]
            scaffold_log.write(f"-- scaffold {i} --\nref1->scaffold:\n{str(a1)}\n")
            scaffold_log.write(f"ref2->scaffold:\n{str(a2)}\n\n")
            scores.append((a1.score, a2.score))

        li = max(range(len(scores)), key=lambda j: scores[j][0])
        ri = max(range(len(scores)), key=lambda j: scores[j][1])
        summary.append(
            {
                "best_left_scaffold": li,
                "best_left_scores": scores[li],
                "best_right_scaffold": ri,
                "best_right_scores": scores[ri],
            }
        )
        print(f"breakpoint {idx}: left={li}({scores[li]}), right={ri}({scores[ri]})")

    scaffold_log.close()

    base = df.drop(columns=["ref1", "ref2", "seq1", "seq2"])
    aug = pd.DataFrame(summary)
    out = pd.concat([base.reset_index(drop=True), aug], axis=1)
    out.to_csv(args.out_table, sep="\t", index=False)
    print(f"Augmented scaffold predictions written to {args.out_table}")

# -------------------------------------------
# Main & CLI
# -------------------------------------------

def main():
    p = argparse.ArgumentParser(description="SV refinement OR scaffold reconstruction")
    p.add_argument("file", help="input TSV")
    p.add_argument(
        "--mode", choices=["refine", "scaffold"], default="refine", help="which pipeline to run"
    )
    p.add_argument(
        "--out-table",
        default="augmented_predictions.tsv",
        help="path for augmented output TSV",
    )
    p.add_argument(
        "--split-log", default="split_read_alignments.txt", help="refine pipeline log"
    )
    p.add_argument(
        "--scaffold-log", default="scaffold_alignments.txt", help="scaffold pipeline log"
    )
    p.add_argument("--outdir", default="out", help="SPAdes output base directory")
    # refine flags
    p.add_argument("-l", "--list", action="store_true", help="list breakpoints")
    p.add_argument(
        "-b", "--breakpoints", nargs="+", type=int, help="which breakpoint indices to process"
    )
    p.add_argument(
        "-v", "--verbose", action="count", default=0, help="verbosity (use -vv for more)"
    )
    # scaffold flags
    p.add_argument("--fasta", help="indexed FASTA for extract_region")

    args = p.parse_args()
    if args.mode == "refine":
        run_refine(args)
    else:
        if not args.fasta:
            p.error("--fasta is required in scaffold mode")
        run_scaffold(args)

if __name__ == "__main__":
    main()
