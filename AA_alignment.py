from Bio import Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import pandas as pd
import subprocess
import re
import argparse
import numpy as np

def check_features(ori, chr1, chr2):
    if chr1 != chr2:
        return "interchromosomal"
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
    if not string:
        return None
    seq = Seq(string)
    rc = seq.reverse_complement()
    return str(rc)

def extract_suffix(cigar):
    # Pattern 1: (num1)S(num2)M
    pattern1 = re.compile(r'(\d+)[M]$')
    pattern2 = re.compile(r'(\d+)[S]$')
    pattern3 = re.compile(r'(\d+)[H]$')

    # Check for first pattern (num1)S(num2)M
    match1 = pattern1.search(cigar)
    if match1:
        return False, match1.group(1)  # Return num1 (the number followed by S or H)
    # Check for first pattern (num1)S(num2)M
    match2 = pattern2.search(cigar)
    if match2:
        return True, match2.group(1)  # Return num1 (the number followed by S or H)
    # Check for first pattern (num1)S(num2)M
    match3 = pattern3.search(cigar)
    if match3:
        return True, match3.group(1)  # Return num1 (the number followed by S or H)
    # If no pattern is matched
    raise ValueError("Cigar string doesn't match (num1)S(num2)M or (num1)M(num2)S\n", cigar)

def format_alignment(arr, query_start, ori, ref_chr, ref_offset=0):
    output = []
    seq_len = 0 # Track sequence length

    for i in range(0, len(arr), 4):  # Process each set of 3 rows
        if i >= len(arr) or arr[i].strip() == "":
            continue

        first_row = arr[i]
        seq = first_row.split()[2]  # Extract the sequence part
        seq_len = len(seq)
        end_number = query_start + seq_len - 1 if ori else query_start - seq_len + 1 # Calculate end number depending on orientation

        first_row = f"ref            {ref_chr}:{query_start} {seq} {end_number}"

        # Align second row
        second_row = arr[i+1]
        query_sequence_start = first_row.index(seq) - 1
        pipes_only = ' ' * query_sequence_start + ''.join([ch if not ch.isdigit() else '' for ch in second_row.strip()])

        # Change first number to 1-indexed
        third_row = arr[i+2]
        ref_start = ref_offset + 1
        ref_seq = third_row.split()[2]  # Extract the sequence part from third row
        ref_len = len(ref_seq)
        ref_end = ref_start + ref_len - 1  # Corrected end number for ref

        ref_row = f"query          {" "*len(ref_chr)} {ref_start} {ref_seq} {ref_end}"

        # Since ref row number is always less than query, calculate extra space
        space_diff = first_row.index(seq) - ref_row.index(ref_seq)
        ref_row = " " * space_diff + ref_row.lstrip()

        output.append(first_row)
        output.append(pipes_only)
        output.append(ref_row)
        output.append('')

        # Update query_start for the next set
        query_start = query_start + seq_len if ori else query_start - seq_len

        # Update the ref_offset for next set
        ref_offset += ref_len

    return "\n".join(output)

def cor_sort(row1, row2):
    chr1 = row1["query_chrom"]
    pos1 = row1["query_pos"]
    chr2 = row2["query_chrom"]
    pos2 = row2["query_pos"]
    chr_num1 = int(chr1[3:])
    chr_num2 = int(chr2[3:])

    if chr_num1 < chr_num2:
        return row1, row2
    elif chr_num1 > chr_num2:
        return row2, row1
    else:
        if pos1 < pos2:
            return row1, row2
        elif pos1 > pos2:
            return row2, row1
        else:
            raise ValueError("Identical read positions!")
        
def rev_cigar(cigar):
    # Use regex to find all the ops of the form 'number+letter'
    ops = re.findall(r'\d+[a-zA-Z]', cigar)
    
    reversed_ops = ops[::-1]
    
    # Join them back into a string
    return ''.join(reversed_ops)

def strand_to_AA(ori, l_ori, r_ori):
    if ori == "++":
        return "+-"
    elif ori == "--":
        return "-+"
    elif ori == "+-":
        return "++"
    elif ori == "-+":
        return "--"
    
    
    # result = None
    # if ori == "++":
    #     result = "+-"
    # elif ori == "--":
    #     result = "-+"
    # elif ori == "+-":
    #     result = "++"
    # elif ori == "-+":
    #     result = "--"
    
    # if l_ori == "-" and r_ori == "+":
    #     result = result.replace('-', '.')
    #     result = result.replace('+', '-')
    #     result = result.replace('.', '+')
    
    # return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize split read junctions and calculate breakpoint details to compare with AA.")
    parser.add_argument("--fa", help="path to hg19 fasta file", type=str, default="./hg19/hg19full.fa")
    parser.add_argument("--aln", help="path to 'alignments.tsv' file from 'alignment.py' output", type=str, default="./alignments.tsv")
    parser.add_argument('--all', action='store_true', help="Show all split reads from at least one AA breakpoint end")
    parser.add_argument('--lrs', action='store_true', help="Show only reads corresponding to low read support AA breakpoints")
    args = parser.parse_args()
    fasta_file = args.fa
    df = pd.read_csv(args.aln, sep="\t")
    # aligner = Align.PairwiseAligner(mode="global", open_gap_score = -10, extend_gap_score = -0.5, match_score = 1.0, mismatch_score = -1.0)
    aligner = Align.PairwiseAligner(mode="local", open_gap_score = -10, extend_gap_score = -5, match_score=2, mismatch_score=-2)
    reg_size = 150

    it = iter(df.iterrows())
    for index, row in it:
        homology = row["AA_homology_seq"]
        if row["split"] == False:
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
        break_hom_len = row["AA_homology_len"]

        query_name = row["query_name"]
        query_is_proper = row["proper_pair"]
        query_read_num = row["read_num"]

        next_row = next(it)[1]

        l_row, r_row = cor_sort(row, next_row)

        l_chr = l_row["query_chrom"]
        l_ori = l_row["query_orientation"]
        l_seq = l_row["query_aln_sub"]
        l_cigar = l_row["query_cigar"]
        l_start = l_row["query_pos"]
        l_end = l_row["query_end"] - 1
        l_len = len(l_seq)

        r_chr = r_row["query_chrom"]
        r_ori = r_row["query_orientation"]
        r_seq = r_row["query_aln_sub"]
        r_cigar = r_row["query_cigar"]
        r_start = r_row["query_pos"]
        r_end = r_row["query_end"] - 1
        r_len = len(r_seq)

        read_row = row if len(row["query_aln_full"]) == 100 else next_row
        read_seq = read_row["query_aln_full"] if read_row["query_orientation"] == "+" else reverse_complement(read_row["query_aln_full"])
        read_ori = read_row["query_orientation"]

        # Check if near both breakpoint ends
        if not args.all:
            if not ((break_chr1 == l_chr and break_chr2 == r_chr) or (break_chr1 == r_chr and break_chr2 == l_chr)):
                continue
            if not (((l_start - reg_size <= break_pos1 <= l_end + reg_size) and (r_start - reg_size <= break_pos2 <= r_end + reg_size)) or
                ((l_start - reg_size <= break_pos2 <= l_end + reg_size) and (r_start - reg_size <= break_pos1 <= r_end + reg_size))):
                continue

        if args.lrs and break_read_support > 3:
            continue
        
        l_aln_seq = l_seq
        l_aln_cigar = l_cigar
        l_aln_start = l_start
        l_aln_end = l_end
        r_aln_seq = r_seq
        r_aln_cigar = r_cigar
        r_aln_start = r_start
        r_aln_end = r_end

        if l_ori == "-":
            l_aln_seq = reverse_complement(l_aln_seq)
            l_aln_cigar = rev_cigar(l_aln_cigar)
            temp = l_aln_start
            l_aln_start = l_aln_end
            l_aln_end = temp
        if r_ori == "-":
            r_aln_seq = reverse_complement(r_aln_seq)
            r_aln_cigar = rev_cigar(r_aln_cigar)
            temp = r_aln_start
            r_aln_start = r_aln_end
            r_aln_end = temp
        #Get orientation and number of clipped bases from both cigars
        hom_len = None
        hom_seq = None

        # Calculate homology length
        # read_seq
        # l_aln_seq
        # r_aln_seq
        # l_len
        # r_len

        l_index = read_seq.index(l_aln_seq)
        r_index = read_seq.index(r_aln_seq)

        seq_list = [(l_index, l_index + l_len), (r_index, r_index + r_len)]
        seq_list.sort(key = lambda x: x[0])
        print(seq_list)

        l_seq_end = seq_list[0][1]
        r_seq_start = seq_list[1][0]
        if l_seq_end > r_seq_start:
            hom_seq = read_seq[r_seq_start:l_seq_end]
            hom_len = len(hom_seq)
        elif l_seq_end == r_seq_start:
            hom_len = 0
        else:
            hom_seq = read_seq[l_seq_end:r_seq_start]
            hom_len = -len(hom_seq)
        
        if read_row['query_orientation'] == "-":
                hom_seq = reverse_complement(hom_seq)
        # Extract reference sequences such that both donor sequences meet in the middle
        region1 = f"{l_chr}:{l_aln_start - (100 - l_len)}-{l_aln_end+10}" if l_ori == "+" else f"{l_chr}:{l_aln_end - 10}-{l_aln_start+ (100 - l_len)}"
        sequence1 = extract_region(fasta_file, region1)
        region2 = f"{r_chr}:{r_aln_start - (100 - r_len)}-{r_aln_end+10}" if r_ori == "+" else f"{r_chr}:{r_aln_end - 10}-{r_aln_start+ (100 - r_len)}"
        sequence2 = extract_region(fasta_file, region2)
        
        # Reverse complement supplementary alignment if necessary to align with primary alignment
        sequence1 = reverse_complement(sequence1) if l_ori == "-" else sequence1
        sequence2 = reverse_complement(sequence2) if r_ori == "-" else sequence2
        full_reference = (sequence1 + sequence2).upper()
        
        print("Fetched split read near AA breakpoint. Displaying breakpoint details:")
        print(row[0:4].T.to_string(header=False))
        print(row[4:10].T.to_string(header=False), "\n")

        print("Split read details:")
        print("Full Reference Sequence:", read_row["query_aln_full"])
        print("Full Query Sequence:", read_seq)
        print(l_row[10:11].T.to_string(header=False))
        print("Left Alignment")
        print(l_row[13:20].T.to_string(header=False).replace("query_pos", "ref_start").replace("query_end", "ref_end  "))
        print("length", 3*"\t", l_len)
        print("Reference Sequence", "\t", l_row["query_aln_sub"])
        print("Query Sequence", 2*"\t", l_aln_seq, "\n")

        print(r_row[10:11].T.to_string(header=False))
        print("Right Alignment")
        print(r_row[13:20].T.to_string(header=False).replace("query_pos", "ref_start").replace("query_end", "ref_end  "))
        print("length", 3*"\t", r_len)
        print("Reference Sequence", "\t", r_row["query_aln_sub"])
        print("Query Sequence", 2*"\t", r_aln_seq, "\n")

        print("Split Read Alignment Visualization:")
        print("Left then Right alignment")
        alignments1 = aligner.align(full_reference, l_aln_seq)
        alignments2 = aligner.align(full_reference, r_aln_seq)

        if read_seq.index(l_aln_seq) != 0:
            alignments1, alignments2 = alignments2, alignments1
            l_chr, r_chr = r_chr, l_chr
            l_aln_start, r_aln_start = r_aln_start, l_aln_start
            l_ori, r_ori = r_ori, l_ori

        for aln in alignments1:
            print(format_alignment(format(aln).split('\n'), l_aln_start, l_ori == "+", l_chr))

        if hom_len != None and hom_len < 0:
            print("Gap detected between donor sequences:", hom_seq)
        
        for aln in alignments2:
            print(format_alignment(format(aln).split('\n'), r_aln_start, r_ori == "+", r_chr, r_seq_start))
        
        strand_ori = l_ori + r_ori
        AA_ori = strand_to_AA(strand_ori, l_ori, r_ori)

        print("Calculated Breakpoint")
        print("chr1", 2*"\t", l_chr)
        print("pos1", 2*"\t", l_start if l_ori == "-" else l_end)
        print("chr2", 2*"\t", r_chr)
        print("pos2", 2*"\t", r_start if r_ori == "+" else r_end)
        print("Strand Orientation", "\t", strand_ori)
        print("\"AA\" Orientation", "\t", AA_ori)
        print("Orientation is based on above order of breakpoint ends!")
        print("features", "\t", check_features(AA_ori, l_chr, r_chr))
        print("homology_length", "\t", hom_len)
        print("homology_sequence", "\t", hom_seq)

        print(80*"-")
        # orientations.append()
