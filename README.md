## Usage
- alignment.py
```
usage: alignment.py [-h] [--sum SUM] [--bam BAM]

Collect reads from the ends of a SV.

options:
  -h, --help  show this help message and exit
  --sum SUM   path to AA amplicon summaries file
  --bam BAM   path to bamfile
```
- AA_alignment.py
```
usage: AA_alignment.py [-h] [--fa FA] [--aln ALN]

Visualize split read junctions and calculate breakpoint details to compare
with AA.

options:
  -h, --help  show this help message and exit
  --fa FA     path to hg19 fasta file
  --aln ALN   path to 'alignments.tsv' file from 'alignment.py' output
```
