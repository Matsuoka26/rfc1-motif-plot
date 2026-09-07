# RFC1 motif plot

This repository contains custom Python scripts used to extract RFC1 repeat-containing long reads and generate repeat motif composition plots.

## Coordinates

Reference genome: hg38

RFC1 repeat core:

```text
chr4:39348425-39348479
```

Surrounding region used for extraction of spanning and soft-clipped reads:

```text
chr4:39348119-39348876
```

## Usage

### 1. Extract RFC1 repeat-containing reads from BAM

```bash
python scripts/extract_rfc1_repeat_reads.py \
  --bam sample.sorted.bam \
  --chrom chr4 \
  --core-start 39348425 \
  --core-end 39348479 \
  --search-start 39348119 \
  --search-end 39348876 \
  --out-prefix sample_RFC1
```

### 2. Generate a repeat motif composition plot

```bash
python scripts/rfc1_motif_plot.py \
  --fastq sample_RFC1.repeat_reads.fastq \
  --out sample_RFC1_motif_plot \
  --motifs AAAGGG,AAAAGG,AAAAG,AAAGG,AAGGG,ACAGG,AAAGC \
  --title "RFC1 repeat motif composition plot" \
  --dpi 400
```

## Notes

The extraction script uses a broader surrounding interval to capture both spanning reads and soft-clipped reads around the RFC1 repeat locus.

The plotting script visualizes sequences as provided in the input FASTQ file. No additional flanking sequence is added during plotting. Sequences that do not match the predefined motifs are shown as OTHER.
