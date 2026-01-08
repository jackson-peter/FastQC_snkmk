# Snakefile for FastQC and MultiQC analysis - Paired-end reads

configfile: "config/config.yaml"

FASTQ_DIR = config["fastq_dir"]
FASTQC_DIR = config.get("fastqc_dir", "fastqc_results")
MULTIQC_DIR = config.get("multiqc_dir", "multiqc_results")
EXTENSION = config.get("extension", "fq.gz")  # Default to fq.gz if not specified

import glob
import os
import re

# Find all files with _1 in them
all_files = glob.glob(f"{FASTQ_DIR}/*_1*.fq.gz") + glob.glob(f"{FASTQ_DIR}/*_1*.fastq.gz")

# Extract sample names by removing _1 and everything after it, then the extension
SAMPLES = []
for f in all_files:
    basename = os.path.basename(f)
    # Match anything before _1, capturing the sample name
    match = re.match(r'^(.+)_1', basename)
    if match:
        sample = match.group(1)
        SAMPLES.append(sample)

# Remove duplicates
SAMPLES = list(set(SAMPLES))

print(f"Found {len(SAMPLES)} samples: {SAMPLES}")

# Rule all - defines the final output
rule all:
    input:
        f"{MULTIQC_DIR}/multiqc_report.html"
        
# Run FastQC on each fastq file
rule fastqc:
    input:
        FASTQ_DIR + "/{sample}." + EXTENSION
    output:
        html = f"{FASTQC_DIR}/{{sample}}_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input}
        """

# Run MultiQC to aggregate all FastQC results
rule multiqc:
    input:
        # Collect all FastQC zip files
        lambda wildcards: [f"{FASTQC_DIR}/{sample}_1P_fastqc.zip" for sample in SAMPLES] + 
                         [f"{FASTQC_DIR}/{sample}_2P_fastqc.zip" for sample in SAMPLES]
    output:
        f"{MULTIQC_DIR}/multiqc_report.html"
    params:
        fastqc_dir = FASTQC_DIR,
        outdir = MULTIQC_DIR
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {params.fastqc_dir} -o {params.outdir}
        """