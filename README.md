# FastQC Quality Control Pipeline

A Snakemake workflow for running FastQC and MultiQC on sequencing data.

## Overview

This pipeline automates quality control analysis of FASTQ files using:
- **FastQC**: Individual quality reports for each FASTQ file
- **MultiQC**: Aggregated report combining all FastQC results

## Features

- Processes all FASTQ files in a specified directory
- Works with any file naming convention (Input file extension is configurable)

## Requirements

- Conda or Mamba
- Snakemake >= 7, < 9

## Installation

1. Clone this repository:
```bash
git clone https://github.com/jackson-peter/Fastqc_snkmk.git
cd Fastqc_snkmk
```

2. Create the Snakemake environment:
```bash
mamba create -n snakemake_qc -c conda-forge -c bioconda "snakemake>=7,<9"
conda activate snakemake_qc
```

## Project Structure

```
.
├── Snakefile              # Main workflow
├── config/
│   └── config.yaml        # Configuration file
├── envs/
│   └── qc.yaml           # Conda environment specification
├── run_snakemake.sh      # Execution script
└── README.md
```

## Configuration

Edit `config/config.yaml` to specify your file paths and extension:

```yaml
fastq_dir: "/path/to/your/fastq/files"
fastqc_dir: "/path/to/output/fastqc"
multiqc_dir: "/path/to/output/multiqc"
extension: "fq.gz"  # or "fastq.gz" for ex
```

## Usage

### Quick Start

```bash
conda activate snakemake_qc
bash run_snakemake.sh
```

### Manual Execution

```bash
# Dry run (see what will be executed)
snakemake --use-conda --conda-frontend mamba --cores 8 -n

# Run the pipeline
snakemake --use-conda --conda-frontend mamba --cores 8 -p
```

### Adjust Resources

Modify the `--cores` parameter based on your available resources:
```bash
snakemake --use-conda --conda-frontend mamba --cores 16 -p
```

## Output

The pipeline generates:

```
fastqc/
├── sample1_fastqc.html
├── sample1_fastqc.zip
├── sample2_fastqc.html
├── sample2_fastqc.zip
└── ...

multiqc/
└── multiqc_report.html  # Combined report for all samples
```

The MultiQC report provides an interactive overview of quality metrics across all your samples.
