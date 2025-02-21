
import pandas as pd
import re
import argparse
import warnings
warnings.filterwarnings("ignore")
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Refine AA SVs.")
    parser.add_argument('-l', "--list", help="shows all breakpoints", action="store_true")
    parser.add_argument('-b', "--breakpoints", type=int, help="list info for these breakpoints (use --list to see breakpoint_number)", nargs="+")
    parser.add_argument('-v', "--verbose", help="show information about specific reads", type=int)
    args = parser.parse_args()
    
    all_reads = pd.read_csv("alignments.tsv", sep="\t")
    no_match_splits = all_reads[all_reads["split"] == True]
    no_match_splits = all_reads[all_reads['query_short'].isin(list(no_match_splits[
        (no_match_splits["break_pos1"] != no_match_splits["query_pos"]) & 
        (no_match_splits["break_pos1"] != no_match_splits["query_end"]) & 
        (no_match_splits["break_pos2"] != no_match_splits["query_pos"]) & 
        (no_match_splits["break_pos2"] != no_match_splits["query_end"])
        ]['query_short']))]
    no_match_splits = no_match_splits[no_match_splits["split"] == True]
    no_match_splits['query_pos'] = no_match_splits.apply(lambda x: x['query_chrom'] + ':' + str(x['query_pos']), axis=1)
    no_match_splits['query_end'] = no_match_splits.apply(lambda x: x['query_chrom'] + ':' + str(x['query_end']), axis=1)


    half_matches = all_reads[(all_reads['split'] == False) & (all_reads['query_cigar'].str.contains('S'))]
    half_matches['read_id'] = half_matches.apply(lambda x: str(x['query_short']) + '_' + str(x['read_num']), axis=1)
    half_matches['query_pos'] = half_matches.apply(lambda x: x['query_chrom'] + ':' + str(x['query_pos']), axis=1)
    half_matches['query_end'] = half_matches.apply(lambda x: x['query_chrom'] + ':' + str(x['query_end']), axis=1)


    break_cands = no_match_splits.groupby("query_short")[["query_short", "query_pos", "query_end"]].agg(set)
    break_cands['all_pos'] = break_cands.apply(lambda x: x['query_pos'].union(x['query_end']), axis=1)
    break_cands = break_cands.drop(['query_short', 'query_pos', 'query_end'], axis=1)


    break_cands['half_matches'] = break_cands.apply(lambda x: half_matches[half_matches["query_pos"].isin(x['all_pos']) | half_matches["query_end"].isin(x['all_pos'])], axis=1)


    all_reads["break_pos1"].value_counts()


    sv = all_reads[all_reads["AA_homology_seq"] == 20950574]


    left = sv[(sv["query_pos"] >= sv["break_pos1"] - 350) & (sv["query_end"] <= sv["break_pos1"] + 350)]
    right = sv[(sv["query_pos"] >= sv["break_pos2"] - 350) & (sv["query_end"] <= sv["break_pos2"] + 350)]


    beg = left["query_pos"].min()
    end = left["query_end"].max()
    left['beg'] = left['query_pos'] - beg


    beg = right["query_pos"].min()
    end = right["query_end"].max()
    right['beg'] = right['query_pos'] - beg


    def refine_step1(reads, certainty = 1):
        beg_reg = r"^\d+(S|H)"
        end_reg = r"\d+(S|H)$"
        filtered_reads = reads[(reads["query_cigar"].str.contains("S")) | (reads["query_cigar"].str.contains("H"))]

        filtered_reads["begin"] = filtered_reads.query_cigar.apply(lambda x: bool(re.search(beg_reg, x)))
        filtered_reads["end"] = filtered_reads.query_cigar.apply(lambda x: bool(re.search(end_reg, x)))

        filtered_reads["sv_end"] = filtered_reads.apply(
        lambda row: row["query_pos"] if row["begin"]
        else row["query_end"] if row["end"]
        else pd.NA, axis=1
        )

        filtered_reads = filtered_reads.dropna(subset=["sv_end"])

        if args.verbose == 2:
            print(filtered_reads.sort_values(by ="sv_end")[['sv_end', 'query_name', 'query_pos', 'query_end', 'query_cigar', 'split', 'query_aln_full']].to_string(), "\n")


        end_cands = filtered_reads["sv_end"].value_counts()
        best = end_cands.index[0] if end_cands.shape[0] > 0 else None

        if args.verbose == 1:
            print(filtered_reads[filtered_reads['sv_end'] == best][['sv_end', 'query_name', 'query_pos', 'query_end', 'query_cigar', 'split', 'query_aln_full']].to_string(), "\n")

        return end_cands, best


    def refine_step2(best):
        contains = lambda read, break_cand: break_cand in range(read["query_pos"] + 1, read["query_end"])
        # for break_cand in cands.index:
        #     display(break_cand, all_reads[all_reads.apply(lambda x: contains(x, break_cand), axis=1)])
        thing = all_reads.loc[all_reads.apply(lambda x: contains(x, best), axis=1)]
        if args.verbose:
            print(thing[['query_name', 'query_pos', 'query_end', 'query_cigar', 'split', 'query_aln_full']].to_string() if len(thing) != 0 else "None", '\n')
        return len(thing)


    pd.set_option('mode.chained_assignment', None)
    svs = all_reads.groupby("break_pos1")
    for group_name, sv in svs:
        left = sv[(sv["query_pos"] >= (sv["break_pos1"] - 350)) & (sv["query_end"] <= (sv["break_pos1"] + 350))]
        right = sv[(sv["query_pos"] >= (sv["break_pos2"] - 350)) & (sv["query_end"] <= (sv["break_pos2"] + 350))]
        if args.list:
            out = svs[["break_chrom1", "break_pos1", "break_chrom2", "break_pos2", "AA_homology_seq"]].first().reset_index(drop=True)
            out.index.name = "breakpoint_index"
            print(out.to_string())
            break
        if args.breakpoints:
            breaks = svs["break_pos1"].first().reset_index(drop=True).iloc[args.breakpoints].to_list()
            if group_name not in breaks:
                continue

        if args.verbose:
            print(f"Reads for {'best' if args.verbose == 1 else ('' if args.verbose == 2 else '')} left SV candidates")
        left_cands, best_left = refine_step1(left)
        if args.verbose:
            print(f"Reads for {'best' if args.verbose == 1 else ('' if args.verbose == 2 else '')} right SV candidates")
        right_cands, best_right = refine_step1(right)

        print("Original AA breakpoint:")
        print(sv.iloc[0][["break_chrom1", "break_pos1", "break_chrom2", "break_pos2", "AA_homology_seq"]].to_string(), '\n')

        print("Best left candidate:", best_left)
        print("Best right candidate:", best_right, '\n')
        
        if args.verbose: 
            print("Overlapping reads with best left candidate:")
        left_ovrs = refine_step2(best_left)
        if args.verbose: 
            print("Overlapping reads with best right candidate:")
        right_ovrs = refine_step2(best_right)

        print("Reads overlapping on best left candidate:", left_ovrs)
        print("Reads overlapping on best right candidate:", right_ovrs)
        print()





