# PanVariants Usage Documentation

## Overview
`PanVariants` is a comprehensive pan-genome variant detection pipeline designed to identify genetic variations across diverse populations. It supports execution on both local machines and Sun Grid Engine (SGE) clusters. Leveraging pan-genome references, the tool enables high-precision detection of Single Nucleotide Variants (SNVs), INDELs (<50bp), Copy Number Variants (CNVs), Short Tandem Repeats (STRs), and other Structural Variants (SVs, ≥50bp), offering robust capabilities for complex genomic analysis.

---

## Quick Start

### Simple Run Command
```bash
python3 PanVariant.py --samplelist samplelist --executor local
```

---
## Install
```bash
git clone https://github.com/MGI-EU/PanVariants.git
```

### Database download
- Here we take the download of databases related to the GRCh38 reference genome as an example; the download process for other reference genomes is similar.
```bash
sh PanVariants_required_files.sh
```

### singularity .sif files download
- Download all .sif files from this URL: https://zenodo.org/records/19848582 and place them in the "sifs" folder.
  
### Required Software
| Software | Version | Installation |
|----------|---------|--------------|
| **Nextflow** | 23.10.0 | Download from: `https://github.com/nextflow-io/nextflow/releases/download/v23.10.0/nextflow-23.10.0-all` <br> Rename to `nextflow` and copy to software directory |
| **Java** | ≥ 17 | Required for Nextflow execution |
| **Singularity** | ≥ 3.8 | Required for containerized execution |
| **Python** | ≥ 3.9.13 | Required for running the Python wrapper script |

### Optional Software
| Software | Description |
|----------|-------------|
| **SGE (Sun Grid Engine)** | Optional cluster management system for distributed computing (required only for `blc` executor) |
---

## Arguments Description

### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--samplelist` | `-s` | Sample information file. Must contain a single line with 3 tab-separated columns:`sample_name` &nbsp; `read1_path` &nbsp; `read2_path` |
| `--executor` | `-ex` | Execution engine options:- `blc`: Run pipeline using a Sun Grid Engine cluster- `local`: Run pipeline using the local machine |

---

### Optional Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--output` | `-o` | `result` | Output directory path |
| `--reference` | `-ref` | `GRCh38` | Reference genome version. Options: `GRCh38` or `CHM13` |
| `--run_cnv` | `-rc` | `yes` | Enable CNV detection (`yes` / `no`) |
| `--run_sv` | `-rs` | `yes` | Enable SV detection (`yes` / `no`) |
| `--run_str` | `-rstr` | `yes` | Enable STR detection (`yes` / `no`) |
| `--mark_dup` | `-md` | `no` | Mark duplicate reads (`yes` / `no`) |
| `--fq_filter` | `-ft` | `no` | Run FastQ filtering (`yes` / `no`) |
| `--split_bam` | `-sb` | `no` | Split BAM files by chromosome and run DeepVariant per chromosome (`yes` / `no`) |
| `--adapter1` | `-ad1` | `AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA` | Adapter 1 sequence |
| `--adapter2` | `-ad2` | `AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG` | Adapter 2 sequence |
| `--re_alignment` | `-ra` | `no` | Perform re-alignment (`yes` / `no`) |
| `--nf_config` | `-c` | - | Path to Nextflow configuration file |
| `--db` | `-d` | - | Database path |
| `--script` | `-sp` | - | Script path |
| `--sif` | `-f` | - | Path to Singularity container image (`.sif`) |
| `--scratch_tmp` | `-scr` | - | Path for Nextflow scratch directive |
| `--queue` | `-q` | - | Queue name for the `qsub -q` parameter. **Note:** Only effective when `--executor` is set to `blc`. |
| `--project` | `-pj` | - | Project name for the `qsub -P` parameter. **Note:** Only effective when `--executor` is set to `blc`. |
| `--ansi_log` | `-an` | `false` | Enable ANSI logging for Nextflow (`true` / `false`) |
| `--help` | `-h` | - | Show help message and exit |

---

## Example Commands

### Run Locally with All Variant Types Enabled
```bash
python3 PanVariant.py \
  --samplelist samples.tsv \
  --executor local \
  --output results \
  --reference GRCh38 \
  --run_cnv yes \
  --run_sv yes \
  --run_str yes
```

### Run on SGE Cluster with Specific Queue and Project
```bash
python3 PanVariant.py \
  --samplelist samples.tsv \
  --executor blc \
  --queue high_mem \
  --project my_project \
  --output cluster_results
```

---

## Notes

- The `--queue` and `--project` parameters are only valid when `--executor` is set to `blc`.
- The sample list file must strictly follow the format: three tab-separated columns with no header row.
- If using Singularity containers, ensure the correct `.sif` file path is provided via `-f`.
- Duplicate marking and FastQ filtering are disabled by default. Set them explicitly to `yes` if required.

---
