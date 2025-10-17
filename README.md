
---

# ecDNA-Alignment

Tools to collect reads around structural-variant (SV) breakpoints and refine AA (AmpliconArchitect) predictions using split-read evidence and optional de-novo scaffolds.

## Prerequisites

* Linux / macOS (bash)
* **git**
* **conda** or **mamba**

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/tshneour/ecDNA-Alignment.git
```
 
### 2. Create conda environment and install requirements
```bash
conda create -n sv-analysis -c conda-forge -c bioconda \
  python>=3.11 spades=4.2.0 samtools pysam biopython pandas numpy natsort

conda activate sv-analysis

```

---

## Quickstart


1. **Setup**

```bash
conda activate sv-analysis
```

2. **Collect reads** for all AA breakpoints into an alignment table + per-breakpoint read pairs:

```
python /path/to/collect.py 350 path/to/AA_summaries/ path/to/sample.bam \
  --strict \
  -v \
  -f alignments.tsv
```

This writes:

* `alignments.tsv` (or `<bam_basename>.tsv` if `-f/--file` not provided)
* gzip’d FASTQs per breakpoint in `fastq/`:

  * `fastq/b_<chr1>_<pos1+1>_<chr2>_<pos2+1>_1.fastq.gz`
  * `fastq/b_<chr1>_<pos1+1>_<chr2>_<pos2+1>_2.fastq.gz`

3. **Run Scaffold and Split Read** and get a combined table:

```bash
python /path/to/refine.py alignments.tsv \
  --mode both \
  --fasta /path/to/genome.fa \
  --out-table final_augmented \
  -v
# Produces final_augmented.tsv
```

> Tip: **Do not** add “.tsv” to `--out-table`; the program appends “.tsv” automatically.

---

## Usage

### `collect.py`

Collect reads from the ends of AA SVs for refinement and write per-breakpoint read pairs (FASTQs) plus a combined TSV of read evidence.

```
usage: collect.py [-h] [-v] [--strict] [-f FILE] refine sum bam
```

**Positional arguments**

* `refine` (int): Radius (bp) around each breakpoint to fetch alignments (e.g., 350).
* `sum` (path): Directory containing AA amplicon summary `.tsv` files.
* `bam` (path): Coordinate-sorted BAM with index (`.bai`).

**Options**

* `-v, --verbose`           Verbose output.
* `--strict`                Only keep read pairs whose alignments fully fall within the region(s) of interest.
* `-f, --file FILE`         Output TSV path (default: `<bam_basename>.tsv`).

**Outputs**

* TSV of reads (`split` and `nonsplit`) with columns including:

  * `break_chrom1`, `break_pos1`, `break_chrom2`, `break_pos2`, `break_sv_type`, `break_read_support`, `break_features`, `break_orientation`, `AA_homology_len`, `AA_homology_seq`, `query_*`, `proper_pair`, `split`, `amplicon`
* Per-breakpoint FASTQs in `fastq/`, gzipped.

### `refine.py`

Refine AA SVs using split reads and/or a de-novo scaffold reconstruction around each breakpoint.

```
usage: refine.py FILE [--mode {split,scaffold,both}]
                     [--out-table OUT] [--split-log PATH]
                     [--scaffold-log PATH] [--outdir DIR]
                     [-l | --list] [-b IDX [IDX ...]] [-v ...]
                     [--fasta FASTA]   # required for scaffold/both
```

**Required**

* `FILE` (path): TSV from `collect.py`.

**Core options**

* `--mode {split,scaffold,both}`   Default: `split`.
* `--out-table OUT`                Basename for output table (program appends `.tsv`).
* `--split-log PATH`               Log of split-read candidate evidence (default: `split_read_alignments.txt`).
* `--scaffold-log PATH`            Log of scaffold alignments per breakpoint (default: `scaffold_alignments.txt`).
* `--outdir DIR`                   SPAdes output base directory (default: `out/`).

**Breakpoint selection & verbosity**

* `-l, --list`                     Print a table of all breakpoints (with indices) and exit.
* `-b, --breakpoints IDX ...`      Only process these `breakpoint_index` values (see `--list`).
* `-v`                             Increase verbosity (`-v`, `-vv`, …).

**Scaffold mode requirement**

* `--fasta FASTA`                  FASTA with `samtools faidx` index present (`.fai`).

**Outputs**

* **Split mode** columns added (per breakpoint):

  * `split_matches`, `sp_left_sv`, `sp_right_sv`,
  * `sp_hom_len` (positive=homology length, negative=insertion length),
  * `hom` (sequence of homology or insertion).
* **Scaffold mode** columns added:

  * `sc_pos1`, `sc_pos2` (estimated refined positions),
  * `sc_hom_len` (Int64; negative=insertion, 0=abutting, positive=homology),
  * `sc_hom` (sequence).
* **Both**: scaffold + split columns combined in one table.

---

## Working assumptions & tips

* **BAM**: Coordinate-sorted, indexed (`.bai` present).
* **AA summaries**: `collect.py` expects `.tsv` files in the summaries folder with AA breakpoints and metadata (e.g., `chrom1/pos1/chrom2/pos2/sv_type/...`). Filenames are used to derive `amplicon` IDs.
* **Mapping quality**: `collect.py` keeps reads with MAPQ > 15 (or mapped status) and writes both split and non-split pairs.
* **FASTQs**: Reconstructed from reads overlapping each pair of breakpoints and written to `fastq/` (gzipped). These are used by scaffold mode.
* **FASTA**: For scaffold mode, index your reference first:

  ```bash
  samtools faidx /path/to/genome.fa
  ```

---
