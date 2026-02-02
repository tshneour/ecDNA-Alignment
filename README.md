# ecDNA-Alignment

Tools to collect reads around structural-variant (SV) breakpoints and refine AA (AmpliconArchitect) predictions using split-read evidence and optional de-novo scaffolds.

## Prerequisites

* Linux / macOS (bash)
* **git**
* **conda** (with the **libmamba** solver enabled; conda ≥ 23 recommended)

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/tshneour/ecDNA-Alignment.git
```

### 2. Create conda environment and install requirements (libmamba)

Use conda with the libmamba solver for faster and more reliable dependency resolution:

```bash
conda create --solver libmamba -n sv-analysis \
  -c conda-forge -c bioconda \
  python>=3.11 spades=4.2.0 samtools pysam biopython pandas numpy natsort

conda activate sv-analysis
```

> Tip: If libmamba is not yet enabled in your conda installation, update conda first:
>
> ```bash
> conda update -n base conda
> ```

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

## Input format

### SV summary (`sum`) files

`collect.py` expects one or more **tab-separated (****`.tsv`****) files** describing structural-variant breakpoints with a fixed set of required columns.

Each TSV in the `sum/` directory represents a collection of SV breakpoints. The filename is used only to derive an `amplicon` identifier (taken as `filename.split("_")[1]`).

#### Required columns

Each input TSV **must** contain the following columns (case-sensitive):

| Column name         | Type   | Description                                                                               |
| ------------------- | ------ | ----------------------------------------------------------------------------------------- |
| `chrom1`            | string | Chromosome of the first breakpoint end (e.g. `chr8`)                                      |
| `pos1`              | int    | 0-based genomic coordinate of the first breakpoint                                        |
| `chrom2`            | string | Chromosome of the second breakpoint end                                                   |
| `pos2`              | int    | 0-based genomic coordinate of the second breakpoint                                       |
| `sv_type`           | string | Structural variant type (e.g. `deletion`, `duplication`, `interchromosomal`, `inversion`) |
| `read_support`      | int    | Number of reads supporting the breakpoint (used for reporting only)                       |
| `features`          | string | Arbitrary annotation or metadata for the breakpoint                                       |
| `orientation`       | string | Breakpoint orientation as a 2-character string (`++`, `--`, `+-`, `-+`)                   |
| `homology_length`   | int    | Length of homology reported for the breakpoint (may be 0)                                 |
| `homology_sequence` | string | Homology or inserted sequence (may be empty)                                              |

Coordinates in the output tables are reported as **1-based**, but the input `pos1` / `pos2` values are treated as **0-based** internally.

---

## Usage

### `collect.py`

Collect reads around SV breakpoints for refinement and write per-breakpoint paired FASTQs plus a combined TSV of read evidence.

```
usage: collect.py [-h] [-v] [--strict] [-f FILE] refine sum bam
```

**Positional arguments**

* `refine` (int): Radius (bp) around each breakpoint to fetch alignments (e.g., 350).
* `sum` (path): Directory containing SV summary `.tsv` files in the format described above.
* `bam` (path): Coordinate-sorted BAM with index (`.bai`).

**Options**

* `-v, --verbose`           Verbose output.
* `--strict`                Only keep read pairs whose alignments fully fall within the region(s) of interest.
* `-f, --file FILE`         Output TSV path (default: `<bam_basename>.tsv`).

**Outputs**

* TSV of reads (`split` and `nonsplit`) with columns including:

  * `break_chrom1`, `break_pos1`, `break_chrom2`, `break_pos2`
  * `break_sv_type`, `break_read_support`, `break_features`, `break_orientation`
  * `query_*` (read-level alignment details)
  * `proper_pair`, `split`, `amplicon`
  * `AA_homology_len`, `AA_homology_seq` (if present in the input SV summary TSVs)

* Per-breakpoint FASTQs in `fastq/`, gzipped:

  * `fastq/b_<chr1>_<pos1>_<chr2>_<pos2>_1.fastq.gz`
  * `fastq/b_<chr1>_<pos1>_<chr2>_<pos2>_2.fastq.gz`

### `refine.py`

Refine SV breakpoints using split-read evidence and/or local de-novo scaffold reconstruction.

```
usage: refine.py FILE [--mode {split,scaffold,both}]
                     [--out-table OUT] [--split-log PATH]
                     [--scaffold-log PATH] [--outdir DIR]
                     [-l | --list] [-b IDX [IDX ...]] [-v ...]
                     [--fasta FASTA]
```

**Required**

* `FILE` (path): TSV produced by `collect.py`.

**Modes**

* `split`: infer refined breakpoint positions and homology using split-read evidence only.
* `scaffold`: perform local assembly (SPAdes) around each breakpoint and align scaffolds to the reference.
* `both`: run scaffold refinement first, then augment with split-read results.

**Notes on AA homology columns**

* `refine.py` does not require `AA_homology_len` / `AA_homology_seq` to be present in the input TSV.

**Outputs**

* **Split mode** columns added:

  * `split_support`, `soft_support`
  * `sp_left_sv`, `sp_right_sv`
  * `sp_hom_len` (positive = homology, negative = insertion)
  * `hom` (homology or insertion sequence)

* **Scaffold mode** columns added:

  * `sc_pos1`, `sc_pos2` (refined breakpoint coordinates)
  * `sc_hom_len` (Int64; negative = insertion, 0 = abutting, positive = homology)
  * `sc_hom` (sequence)

* **Both**: combined scaffold and split-read columns in a single TSV.

---

## Working assumptions & tips

* **BAM**: Coordinate-sorted, indexed (`.bai` present).
* **SV summaries**: Input TSVs must conform exactly to the column specification above.
* **Mapping quality**: Reads with MAPQ > 15 (or mapped status) are retained.
* **FASTQs**: Reconstructed from reads overlapping each breakpoint pair and used by scaffold mode.
* **FASTA**: Required for scaffold mode; must be indexed first:

  ```bash
  samtools faidx /path/to/genome.fa
  ```
