import pandas as pd
import re
import argparse
import warnings

warnings.filterwarnings("ignore")
print_columns = [
    "query_name",
    "proper_pair",
    "query_pos",
    "query_end",
    "query_cigar",
    "split",
    "query_aln_full",
]
print_columns2 = [
    "query_short",
    "query_chrom",
    "query_pos",
    "query_cigar",
    "does_clip_match",
    "rev_clipped",
    "query_aln_sub",
]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Refine AA SVs.")
    parser.add_argument("file", type=str, help="alignments file to refine")
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

    # ------------- Define Functions ------------------

    def refine_step1(reads):
        reads = reads[reads["query_cigar"].str.contains(r"[SH]")].copy()
        reads["begin"] = reads["query_cigar"].str.contains(r"^\d+[SH]", regex=True)
        reads["end"] = reads["query_cigar"].str.contains(r"\d+[SH]$", regex=True)

        reads["sv_end"] = reads["query_pos"].where(reads["begin"])
        reads["sv_end"].fillna(reads["query_end"].where(reads["end"]), inplace=True)

        # Drop rows where 'sv_end' is still NaN
        reads.dropna(subset=["sv_end"], inplace=True)

        if args.verbose == 2:
            print(
                reads.sort_values(by="sv_end")[["sv_end", *print_columns]].to_string(
                    index=False
                ),
                "\n",
            )

        end_cands = reads["sv_end"].value_counts()
        best = end_cands.index[0] if end_cands.shape[0] > 0 else None

        if args.verbose == 1:
            verb_1 = reads[reads["sv_end"] == best][["sv_end", *print_columns]]
            disc = verb_1[verb_1["proper_pair"] == "Discordant"]
            mates = all_reads[
                all_reads["query_name"].isin(disc.query_name)
                & ~all_reads["query_aln_full"].isin(disc.query_aln_full)
            ][print_columns]
            print(verb_1.to_string(index=False), "\n")
            print()
            print("Mates of Discordant Reads from above Table: ")
            print(mates.to_string(index=False), "\n")

        return end_cands, best, reads

    def refine_step2(best):
        def contains(read, break_cand):
            return break_cand in range(read["query_pos"] + 1, read["query_end"])

        thing = all_reads.loc[all_reads.apply(lambda x: contains(x, best), axis=1)]
        if args.verbose:
            print(
                thing[print_columns].to_string(index=False) if thing else "None",
                "\n",
            )
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

        left["clip_len"] = left.apply(
            lambda x: len(x["query_aln_full"]) - len(x["query_aln_sub"]), axis=1
        )
        left["clipped"] = left.apply(
            lambda x: (
                x["query_aln_full"][0 : x["clip_len"]]
                if x["begin"]
                else x["query_aln_full"][-x["clip_len"] or len(x["query_aln_full"]) :]
            ),
            axis=1,
        )
        right["clip_len"] = right.apply(
            lambda x: len(x["query_aln_full"]) - len(x["query_aln_sub"]), axis=1
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
                left_df["does_clip_match"] = left_df["rev_clipped"].apply(
                    lambda x: last_nohomo.startswith(x)
                )
                left_df.loc[left_df["split"], "does_clip_match"] = left_df.loc[
                    left_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (right_df.loc[right_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["does_clip_match"]
                    ),
                    axis=1,
                )
                right_df["rev_clipped"] = right_df["clipped"].where(
                    right_df["begin"], rev_comp(right_df["clipped"])
                )
                right_df["does_clip_match"] = right_df["rev_clipped"].apply(
                    lambda x: first_nohomo.endswith(x)
                )
                right_df.loc[right_df["split"], "does_clip_match"] = right_df.loc[
                    right_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (left_df.loc[left_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["does_clip_match"]
                    ),
                    axis=1,
                )

                hom_sum_left = left_df["does_clip_match"].sum()
                hom_sum_right = right_df["does_clip_match"].sum()

                # Check for insertion
                first_row = left_df.loc[left_df["rev_clipped"].str.len().idxmax()]
                last_row = right_df.loc[right_df["rev_clipped"].str.len().idxmax()]
                first_clipped = first_row["rev_clipped"]
                last_clipped = last_row["rev_clipped"]
                insertion = get_homology(last_clipped, first_clipped)
                insertion_len = len(insertion)
                left_df["does_clip_match"] = (
                    left_df["rev_clipped"]
                    .str.slice(start=insertion_len)
                    .apply(lambda x: last.startswith(x))
                )
                left_df.loc[left_df["split"], "does_clip_match"] = left_df.loc[
                    left_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (right_df.loc[right_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["does_clip_match"]
                    ),
                    axis=1,
                )
                right_df["does_clip_match"] = (
                    right_df["rev_clipped"]
                    .str.slice(stop=-insertion_len)
                    .apply(lambda x: first.endswith(x))
                )
                right_df.loc[right_df["split"], "does_clip_match"] = right_df.loc[
                    right_df["split"]
                ].apply(
                    lambda x: (
                        x["query_short"]
                        in (left_df.loc[left_df["split"], "query_short"]).values
                        if "H" in x["query_cigar"]
                        else x["does_clip_match"]
                    ),
                    axis=1,
                )
                ins_sum_left = left_df["does_clip_match"].sum()
                ins_sum_right = right_df["does_clip_match"].sum()

                if hom_sum_left == len(left_df) and hom_sum_right == len(right_df):
                    print("Yes! All soft clips match!")
                    if best_left is not None:
                        print(
                            "Hmmmm.... Multiple candidates with matching soft clips",
                            f"Namely left: {left_group} right: {right_group} homology: {homology}",
                        )
                    best_left = left_group
                    best_right = right_group
                    best_homology = homology if first_row["end"] else rev_comp(homology)
                elif hom_sum_left > 0 or hom_sum_right > 0:
                    print(
                        "Hmmmmm...... only some soft clips matched",
                        f"Namely left: {hom_sum_left}/{len(left_df)} right: {hom_sum_right}/{len(right_df)}",
                    )
                    print(
                        "More info:",
                        f"Namely left: {left_group} right: {right_group} homology: {homology}",
                    )
                    print(left_df[print_columns2].to_string(index=False))
                    print(right_df[print_columns2].to_string(index=False))
                    print()

                if insertion_len == 0:
                    continue
                if ins_sum_left == len(left_df) and ins_sum_right == len(right_df):
                    print()
                    print("Trying as insertion")
                    print("Yes! All soft clips match!")
                    if best_left is not None:
                        print(
                            "Hmmmm.... Multiple candidates with matching soft clips",
                            f"Namely left: {left_group} right: {right_group} insertion: {insertion}",
                        )
                    best_left_ins = left_group
                    best_right_ins = right_group
                    best_insertion = (
                        insertion if first_row["end"] else rev_comp(insertion)
                    )
                elif ins_sum_left > 0 or ins_sum_right > 0:
                    print()
                    print("Trying as insertion")
                    print(
                        "Hmmmmm...... only some soft clips matched",
                        f"Namely left: {ins_sum_left}/{len(left_df)} right: {ins_sum_right}/{len(right_df)}",
                    )
                    print(
                        "More info:",
                        f"Namely left: {left_group} right: {right_group} insertion: {insertion}",
                    )
                    print(left_df[print_columns2].to_string(index=False))
                    print(right_df[print_columns2].to_string(index=False))
                    print()

        return (
            best_left,
            best_right,
            best_homology,
            best_left_ins,
            best_right_ins,
            best_insertion,
        )

    def rev_comp(series):
        trans = str.maketrans("ACGT", "TGCA")
        return (
            series.str[::-1].str.translate(trans)
            if isinstance(series, pd.Series)
            else series[::-1].translate(trans)
        )

    def get_homology(first, last):
        for i in range(len(first), -1, -1):
            if last.startswith(first[-i:]):
                return first[-i:]
        return ""

    # ------------- Do refinement ------------------

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
        left_cands, best_left, left_groups = refine_step1(left)
        if args.verbose:
            print(
                f"Reads for {'best' if args.verbose == 1 else ('' if args.verbose == 2 else '')} right SV candidates"
            )
        right_cands, best_right, right_groups = refine_step1(right)

        (
            new_best_left,
            new_best_right,
            homology,
            new_best_insL,
            new_best_insR,
            insertion,
        ) = check_overlap(left_groups, right_groups)
        print(
            f"New refinement (Homology) -> left: {new_best_left} right: {new_best_right} homology: {homology}"
        )
        print(
            f"New refinement (Insertion) -> left: {new_best_insL} right: {new_best_insR} insertion: {insertion}"
        )

        print("Original AA breakpoint:")
        print(
            sv.iloc[0][
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
