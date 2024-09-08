from Bio import Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import pandas as pd
import subprocess
import re
import argparse

def check_features(ori):
    if ori == "++":
        return "foldback"
    elif ori == "--":
        return "inversion"
    elif ori == "+-":
        return "deletion-like"
    elif ori == "-+":
        return "duplication-like"
    else:
        return "Unknown"

def extract_region(fasta_file, region):
    """
    Extracts a specific region from a FASTA file using samtools.

    :param fasta_file: Path to the indexed FASTA file.
    :param region: Region in the format "chr:start-end" (e.g., "chr18:180-280").
    :return: The sequence as a string.
    """
    result = subprocess.run(
        ["samtools", "faidx", fasta_file, region],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if result.returncode == 0:
        # The sequence will be in the second line, as the first line contains the header
        return ''.join(result.stdout.splitlines()[1:])
    else:
        raise Exception(f"Error: {result.stderr.strip()}")

def reverse_complement(string):
    seq = Seq(string)
    rc = seq.reverse_complement()
    return str(rc)

def extract_affix(cigar):
    # Pattern 1: (num1)S(num2)M
    pattern1 = re.compile(r'^(\d+)[SH]')
    # Pattern 2: (num1)M(num2)S
    pattern2 = re.compile(r'(\d+)[SH]$')

    # Check for first pattern (num1)S(num2)M
    match1 = pattern1.search(cigar)
    if match1:
        return False, match1.group(1)  # Return num1 (the number followed by S)
    
    # Check for second pattern (num1)M(num2)S
    match2 = pattern2.search(cigar)
    if match2:
        return True, match2.group(1)  # Return num2 (the number followed by S)
    
    # If no pattern is matched
    raise ValueError("Cigar string doesn't match (num1)S(num2)M or (num1)M(num2)S\n", cigar)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize split read junctions and calculate breakpoint details to compare with AA.")
    parser.add_argument("--fa", help="path to hg19 fasta file", type=str, default="./hg19/hg19full.fa")
    parser.add_argument("--aln", help="path to 'alignments.tsv' file from 'alignment.py' output", type=str, default="./alignments.tsv")
    args = parser.parse_args()
    fasta_file = args.fa
    df = pd.read_csv(args.aln, sep="\t")
    # aligner = Align.PairwiseAligner(mode="global", open_gap_score = -10, extend_gap_score = -0.5, match_score = 1.0, mismatch_score = -1.0)
    aligner = Align.PairwiseAligner(mode="local", open_gap_score = -10, extend_gap_score = -5, match_score=2, mismatch_score=-2)

    it = iter(df.iterrows())
    for index, row in it:
        homology = row["break_homology_seq"]
        if row["split"] == False or pd.isna(homology):
            continue
        
        # Collect AA breakpoint details
        break_chr1 = row["break_chrom1"]
        break_pos1 = row["break_pos1"]
        break_chr2 = row["break_chrom2"]
        break_pos2 = row["break_pos2"]
        break_ori = row["break_orientation"]
        break_sv_type = row["break_sv_type"]
        break_read_support = row["break_read_support"]
        break_features = row["break_features"]
        break_hom_len = row["break_homology_len"]

        query_name = row["query_name"]
        query_is_proper = row["proper_pair"]
        query_read_num = row["read_num"]

        second_row = next(it)[1]

        if 'H' in row["query_cigar"]:
            temp = row
            row = second_row
            second_row = temp

        first_chr = row["query_chrom"]
        first_ori = row["query_orientation"]
        # first_aln = Seq(row["query_aln_sub"]) if row["orientation"][0]=="+" else Seq(row["query_aln_sub"]).reverse_complement()
        first_aln = row["query_aln_sub"]
        first_cigar = row["query_cigar"]
        first_pos = row["query_pos"]
        first_end = row["query_end"]
        first_len = len(first_aln)

        second_chr = second_row["query_chrom"]
        second_ori = second_row["query_orientation"]
        # second_aln = Seq(second_row["query_aln_sub"]) if row["orientation"][1]=="+" else Seq(second_row["query_aln_sub"]).reverse_complement()
        second_aln = second_row["query_aln_sub"]
        second_cigar = second_row["query_cigar"]
        second_pos = second_row["query_pos"]
        second_end = second_row["query_end"]
        second_len = len(second_aln)

        primary_aln = row["query_aln_full"]
        primary_len = len(primary_aln)

        #Get orientation and number of clipped bases from both cigars
        pr_ori, pr_affix = extract_affix(first_cigar)
        sup_ori, sup_affix = extract_affix(second_cigar)
        pr_affix = int(pr_affix)
        sup_affix = int(sup_affix)
        hom_len = None
        hom_seq = None

        # Calculate homology length
        if pr_ori:
                sup_aln = reverse_complement(second_aln) if sup_ori else second_aln
                sup_aln_len = len(sup_aln)
                suffix_len = pr_affix
                hom_len = (primary_len - suffix_len) - (primary_len - sup_aln_len)
                if hom_len > 0:
                    hom_seq = primary_aln[primary_len - sup_aln_len: primary_len - suffix_len]
                elif hom_len < 0:
                    hom_seq = primary_aln[primary_len - suffix_len: primary_len - sup_aln_len]
        else:
            sup_aln = second_aln if sup_ori else reverse_complement(second_aln)
            sup_aln_len = len(sup_aln)
            prefix_len = pr_affix
            hom_len = sup_aln_len - prefix_len
            if hom_len > 0:
                hom_seq = primary_aln[prefix_len: sup_aln_len]
            elif hom_len < 0:
                hom_seq = primary_aln[sup_aln_len: prefix_len]

        # Extract reference sequences such that both donor sequences meet in the middle
        region1 = f"{first_chr}:{first_pos - (100 - first_len)}-{first_pos+first_len+9}" if pr_ori else f"{first_chr}:{first_pos - 10}-{first_pos+99}"
        sequence1 = extract_region(fasta_file, region1)
        region2 = f"{second_chr}:{second_pos-10}-{second_pos+99}" if ((not pr_ori and not sup_ori) or (pr_ori and not sup_ori)) else f"{second_chr}:{second_pos - (100 - second_len)}-{second_pos+second_len+9}"
        sequence2 = extract_region(fasta_file, region2)
        
        # Reverse complement supplementary alignment if necessary to align with primary alignment
        second_aln = reverse_complement(second_aln) if (pr_ori and sup_ori) or (not pr_ori and not sup_ori) else second_aln
        sequence2 = reverse_complement(sequence2) if (pr_ori and sup_ori) or (not pr_ori and not sup_ori) else sequence2
        full_reference = (sequence1 + sequence2).upper() if pr_ori else (sequence2 + sequence1).upper()
        
        print("Fetched split read near AA breakpoint. Displaying breakpoint details:")
        print(row[0:4].T.to_string(header=False))
        print(row[4:10].T.to_string(header=False), "\n")

        print("Split read details:")
        print(row[10:11].T.to_string(header=False))
        print("Primary Alignment" if 'S' in first_cigar else "Supplementary Alignment")
        print(row[13:20].T.to_string(header=False))
        print("length", 3*"\t", first_len)
        print("Raw Sequence", 2*"\t", row["query_aln_sub"], "\n")
        print(second_row[10:11].T.to_string(header=False))
        print("Primary Alignment" if 'S' in second_cigar else "Supplementary Alignment")
        print(second_row[13:20].T.to_string(header=False))
        print("length", 3*"\t", second_len)
        print("Raw Sequence", 2*"\t", second_row["query_aln_sub"], "\n")

        print("Split Read Alignment Visualization:")
        print("Primary first then Sup" if pr_ori else "Sup first then Primary")
        alignments1 = aligner.align(full_reference, first_aln if pr_ori else second_aln)
        alignments2 = aligner.align(full_reference, first_aln if not pr_ori else second_aln)
        for aln in alignments1:
            print(aln)
        if hom_len < 0:
            print("Gap detected between donor sequences:", hom_seq)
        for aln in alignments2:
            print(aln)
        
        calc_ori = ("+" if pr_ori else "-") + ("+" if sup_ori else "-")
        print("Calculated Breakpoint")
        print("chr1", 2*"\t", first_chr)
        print("pos1", 2*"\t", first_pos if calc_ori[0] == "-" else first_end)
        print("chr2", 2*"\t", second_chr)
        print("pos2", 2*"\t", second_pos if calc_ori[1] == "-" else second_end)
        print("orientation", "\t", calc_ori)
        print("Orientation is based on above order of breakpoint ends!")
        print("features", "\t", check_features(calc_ori))
        print("homology_length", "\t", hom_len)
        print("homology_sequence", "\t", hom_seq)

        print("\n")
