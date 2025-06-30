#!/usr/bin/env python3
import argparse
import warnings
import pandas as pd

# common
warnings.filterwarnings("ignore")
pd.set_option("mode.chained_assignment", None)

# -------------------------------------------
# Refine pipeline (your original code)
# -------------------------------------------

print_columns = [
    "query_name", "proper_pair", "query_pos", "query_end",
    "query_cigar", "split", "query_aln_full",
]
print_columns2 = [
    "query_short", "query_chrom", "query_pos", "query_cigar",
    "hom_clip_match", "ins_clip_match", "rev_clipped", "query_aln_sub",
]

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

def refine_step1(reads, all_reads, verbose):
    reads = reads[reads["query_cigar"].str.contains(r"[SH]")].copy()
    reads["begin"] = reads["query_cigar"].str.contains(r"^\d+[SH]", regex=True)
    reads["end"] = reads["query_cigar"].str.contains(r"\d+[SH]$", regex=True)

    reads["sv_end"] = reads["query_pos"].where(reads["begin"])
    reads["sv_end"].fillna(reads["query_end"].where(reads["end"]), inplace=True)
    reads.dropna(subset=["sv_end"], inplace=True)

    if verbose == 2:
        print(
            reads.sort_values(by="sv_end")[["sv_end", *print_columns]]
            .to_string(index=False), "\n"
        )

    end_cands = reads["sv_end"].value_counts()
    best = end_cands.index[0] if len(end_cands) > 0 else None

    if verbose == 1 and best is not None:
        verb_1 = reads[reads["sv_end"] == best][["sv_end", *print_columns]]
        disc = verb_1[verb_1["proper_pair"] == "Discordant"]
        mates = all_reads[
            all_reads["query_name"].isin(disc.query_name)
            & ~all_reads["query_aln_full"].isin(disc.query_aln_full)
        ][print_columns]
        print(verb_1.to_string(index=False), "\n")
        print("Mates of Discordant Reads:")
        print(mates.to_string(index=False), "\n")

    return best, reads

def refine_step2(best, all_reads, verbose):
    def contains(read):
        return best in range(read["query_pos"] + 1, read["query_end"])
    thing = all_reads.loc[all_reads.apply(contains, axis=1)]
    if verbose:
        print(thing[print_columns].to_string(index=False) if len(thing) else "None", "\n")
    return len(thing)

def check_overlap(left, right):

        # Eliminate unhelpful split reads
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

        left["clip_len"] = (
            left["query_aln_full"].str.len() - left["query_aln_sub"].str.len()
        )
        left["clipped"] = left.apply(
            lambda x: (
                x["query_aln_full"][0 : x["clip_len"]]
                if x["begin"]
                else x["query_aln_full"][-x["clip_len"] or len(x["query_aln_full"]) :]
            ),
            axis=1,
        )
        right["clip_len"] = (
            right["query_aln_full"].str.len() - right["query_aln_sub"].str.len()
        )
        right["clipped"] = right.apply(
            lambda x: (
                x["query_aln_full"][0 : x["clip_len"]]
                if x["begin"]
                else x["query_aln_full"][-x["clip_len"] or len(x["query_aln_full"]) :]
            ),
            axis=1,
        )
        left_groups = left.groupby("sv_end")
        right_groups = right.groupby("sv_end")

        results = []

        for left_group, left_df in left_groups:
            for right_group, right_df in right_groups:
                first_row = left_df.loc[left_df["query_aln_sub"].str.len().idxmax()]
                last_row = right_df.loc[right_df["query_aln_sub"].str.len().idxmax()]
                first = (
                    first_row["query_aln_sub"]
                    if first_row["end"]
                    else rev_comp(first_row["query_aln_sub"])
                )
                last = (
                    last_row["query_aln_sub"]
                    if last_row["begin"]
                    else rev_comp(last_row["query_aln_sub"])
                )

                # Check for homology
                homology = get_homology(first, last)
                hom_len = len(homology)
                first_nohomo = first[:-hom_len]
                last_nohomo = last[hom_len:]
                left_df["rev_clipped"] = left_df["clipped"].where(
                    left_df["end"], rev_comp(left_df["clipped"])
                )
                left_df["hom_clip_match"] = left_df["rev_clipped"].apply(
                    lambda x: last_nohomo.startswith(x)
                )
                left_df.loc[left_df["split"], "hom_clip_match"] = left_df.loc[
                    left_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (right_df.loc[right_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["hom_clip_match"]
                    ),
                    axis=1,
                )
                right_df["rev_clipped"] = right_df["clipped"].where(
                    right_df["begin"], rev_comp(right_df["clipped"])
                )
                right_df["hom_clip_match"] = right_df["rev_clipped"].apply(
                    lambda x: first_nohomo.endswith(x)
                )
                right_df.loc[right_df["split"], "hom_clip_match"] = right_df.loc[
                    right_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (left_df.loc[left_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["hom_clip_match"]
                    ),
                    axis=1,
                )

                hom_sum_left = left_df["hom_clip_match"].sum()
                hom_sum_right = right_df["hom_clip_match"].sum()

                # Check for insertion
                first_row = left_df.loc[left_df["rev_clipped"].str.len().idxmax()]
                last_row = right_df.loc[right_df["rev_clipped"].str.len().idxmax()]
                first_clipped = first_row["rev_clipped"]
                last_clipped = last_row["rev_clipped"]
                insertion = get_homology(last_clipped, first_clipped)
                insertion_len = len(insertion)
                left_df["ins_clip_match"] = (
                    left_df["rev_clipped"]
                    .str.slice(start=insertion_len)
                    .apply(lambda x: last.startswith(x))
                )
                left_df.loc[left_df["split"], "ins_clip_match"] = left_df.loc[
                    left_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (right_df.loc[right_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["ins_clip_match"]
                    ),
                    axis=1,
                )
                right_df["ins_clip_match"] = (
                    right_df["rev_clipped"]
                    .str.slice(stop=-insertion_len)
                    .apply(lambda x: first.endswith(x))
                )
                right_df.loc[right_df["split"], "ins_clip_match"] = right_df.loc[
                    right_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (left_df.loc[left_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["ins_clip_match"]
                    ),
                    axis=1,
                )
                ins_sum_left = left_df["ins_clip_match"].sum()
                ins_sum_right = right_df["ins_clip_match"].sum()

                split_matches = (
                    left_df["split"].sum() +
                    right_df["split"].sum()
                )

                results.append(
                    pd.DataFrame(
                        {
                            "left_sv": left_group,
                            "right_sv": right_group,
                            "hom_sum_left": hom_sum_left,
                            "hom_sum_right": hom_sum_right,
                            "homology": homology,
                            "ins_sum_left": ins_sum_left,
                            "ins_sum_right": ins_sum_right,
                            "insertion": insertion,
                            "left_len": len(left_df),
                            "right_len": len(right_df),
                            "split_matches": split_matches,
                            "left": (left_df.copy(deep=True),),
                            "right": (right_df.copy(deep=True),),
                        },
                        index=[0],
                    )
                )

        results = pd.concat(results)
        results["hom_%"] = (
            (results["hom_sum_left"] + results["hom_sum_right"])
            / (results["left_len"] + results["right_len"])
        )
        results["ins_%"] = (
            (results["ins_sum_left"] + results["ins_sum_right"])
            / (results["left_len"] + results["right_len"])
        )
        results["total_%"] = results[["hom_%", "ins_%"]].max(axis=1)
        results["split_matches"] = results["split_matches"]
        results["total_reads"] = results["left_len"] + results["right_len"]
        results = results.sort_values(by=["total_%", "split_matches", "total_reads"], ascending=False)

        return results.head(3)[["left_sv", "right_sv", "hom_%", "ins_%", "total_%", "total_reads", "split_matches", "homology", "insertion", "left", "right"]]


def run_refine(args):
    all_reads = pd.read_csv(args.file, sep="\t")
    svs = all_reads.groupby("break_pos1")

    if args.list:
        out = (
            svs[["break_chrom1","break_pos1","break_chrom2","break_pos2","AA_homology_seq"]]
            .first().reset_index(drop=True)
        )
        out.index.name = "breakpoint_index"
        print(out.to_string())
        return

    breaks = None
    if args.breakpoints:
        idxs = svs["break_pos1"].first().reset_index(drop=True)
        breaks = set(idxs.iloc[args.breakpoints].to_list())

    for bp, sv in svs:
        if breaks and bp not in breaks:
            continue

        left = sv[
            (sv["query_pos"] >= bp - 350)
            & (sv["query_end"] <= bp + 350)
        ]
        right = sv[
            (sv["query_pos"] >= sv["break_pos2"].iloc[0] - 350)
            & (sv["query_end"] <= sv["break_pos2"].iloc[0] + 350)
        ]

        best_left, left_grp = refine_step1(left, all_reads, args.verbose)
        best_right, right_grp = refine_step1(right, all_reads, args.verbose)

        results = check_overlap(left_grp, right_grp)
        for i in range(len(results)-1, -1, -1):
            row = results.iloc[i]
            print("_" * 100)
            print(row.drop(["left","right"]).to_string())
            print("\nLeft reads:")
            print(row["left"][print_columns2].to_string(index=False))
            print("\nRight reads:")
            print(row["right"][print_columns2].to_string(index=False))
            print("_" * 100, "\n\n")

        print("Original AA breakpoint:")
        bp0 = sv.iloc[0]
        print(
            bp0[["break_chrom1","break_pos1","break_chrom2","break_pos2",
                 "AA_homology_seq","break_orientation"]]
            .to_string(), "\n"
        )

# -------------------------------------------
# Scaffold pipeline (from your notebook)
# -------------------------------------------
from Bio import Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import subprocess, shutil

aligner = Align.PairwiseAligner(
    mode="local", open_gap_score=-10, extend_gap_score=-5,
    match_score=2, mismatch_score=-2
)

def extract_region(fasta_file, region):
    res = subprocess.run(
        ["samtools","faidx",fasta_file,region],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if res.returncode:
        raise RuntimeError(res.stderr)
    return "".join(res.stdout.splitlines()[1:])

def generate_scaffolds(fq1, fq2, out_dir="out"):
    shutil.rmtree(out_dir, ignore_errors=True)
    cmd = [
        "python","./SPAdes-4.2.0-Linux/bin/spades.py","--meta",
        "--pe1-1", fq1, "--pe1-2", fq2, "-o", out_dir
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode:
        print(res.stdout)
        print(res.stderr)
        raise RuntimeError(res.stderr)
    seqs = []
    with open(f"{out_dir}/scaffolds.fasta") as F:
        for chunk in F.read().split(">")[1:]:
            seqs.append("".join(chunk.splitlines()[1:]))
    return seqs

def run_scaffold(args):
    # load breakpoints
    df = pd.read_csv(args.file, sep="\t")
    # build ref windows
    df["ref1"] = df.break_chrom1 + ":" + (df.break_pos1-350).astype(str) + "-" + (df.break_pos1+350).astype(str)
    df["ref2"] = df.break_chrom2 + ":" + (df.break_pos2-350).astype(str) + "-" + (df.break_pos2+350).astype(str)
    df["seq1"] = df.ref1.apply(lambda r: extract_region(args.fasta, r))
    df["seq2"] = df.ref2.apply(lambda r: extract_region(args.fasta, r))

    # assume fastqs in ./fastq/, names b_chr1_pos1_chr2_pos2_1.fastq.gz & _2
    def fq_names(row):
        base = f"b_{row.break_chrom1}:{row.break_pos1}_{row.break_chrom2}:{row.break_pos2}"
        return f"./fastq/{base}_1.fastq.gz", f"./fastq/{base}_2.fastq.gz"

    scaffs = []
    for _, row in df.iterrows():
        fq1, fq2 = fq_names(row)
        scaffs.append(generate_scaffolds(fq1, fq2, out_dir=args.outdir))
    df["scaffolds"] = scaffs

    # align & score
    results = []
    for idx, row in df.iterrows():
        scores = [(aligner.align(row.seq1.upper(), sc)[0].score,
                   aligner.align(row.seq2.upper(), sc)[0].score)
                  for sc in row.scaffolds]
        best_left = max(range(len(scores)), key=lambda i: scores[i][0])
        best_right= max(range(len(scores)), key=lambda i: scores[i][1])
        print(f"breakpoint {idx}: best_left={best_left} (score {scores[best_left]}), "
              f"best_right={best_right} (score {scores[best_right]})")

# -------------------------------------------
# Main & CLI
# -------------------------------------------
def main():
    p = argparse.ArgumentParser(description="AA SV refinement OR scaffold reconstruction")
    p.add_argument("file", help="input TSV")
    sub = p.add_argument_group("common options")
    sub.add_argument("--mode", choices=["refine","scaffold"], default="refine",
                     help="which pipeline to run")
    # refine flags
    sub.add_argument("-l","--list", action="store_true", help="list breakpoints")
    sub.add_argument("-b","--breakpoints", type=int, nargs="+",
                     help="which breakpoint indices to do")
    sub.add_argument("-v","--verbose", action="count", default=0,
                     help="verbosity (use -vv for full)")
    # scaffold flags
    sub.add_argument("--fasta", help="indexed FASTA for extract_region")
    sub.add_argument("--outdir", default="out", help="SPAdes output directory")

    args = p.parse_args()
    if args.mode == "refine":
        run_refine(args)
    else:
        if not args.fasta:
            p.error("--fasta is required in scaffold mode")
        run_scaffold(args)

if __name__ == "__main__":
    main()
