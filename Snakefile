# Snakefile for FastQC and MultiQC analysis - Paired-end reads

# Configuration
FASTQ_DIR = "fastq"  # Directory containing your .fq.gz or .fastq.gz files
FASTQC_DIR = "fastqc_results"
MULTIQC_DIR = "multiqc_results"

# Get all fastq.gz files and extract sample names
# Assumes naming convention: sample_R1.fq.gz and sample_R2.fq.gz
# Or: sample_1.fq.gz and sample_2.fq.gz
import glob
import os
import re

# Find all R1/1 files
r1_files = glob.glob(f"{FASTQ_DIR}/*_R1.fq.gz") + glob.glob(f"{FASTQ_DIR}/*_R1.fastq.gz") + \
           glob.glob(f"{FASTQ_DIR}/*_1.fq.gz") + glob.glob(f"{FASTQ_DIR}/*_1.fastq.gz")

# Extract sample names (remove _R1/_1 and extension)
SAMPLES = []
for f in r1_files:
    basename = os.path.basename(f)
    sample = re.sub(r'_(R)?1\.(fq|fastq)\.gz$', '', basename)
    SAMPLES.append(sample)

# Rule all - defines the final output
rule all:
    input:
        f"{MULTIQC_DIR}/multiqc_report.html"

# Run FastQC on R1 files
rule fastqc_r1:
    input:
        r1 = FASTQ_DIR + "/{sample}_R1.fq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_R1_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_R1_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1}
        """

# Run FastQC on R2 files
rule fastqc_r2:
    input:
        r2 = FASTQ_DIR + "/{sample}_R2.fq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_R2_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_R2_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r2}
        """

# Alternative rules for _1/_2 naming convention
rule fastqc_1:
    input:
        r1 = FASTQ_DIR + "/{sample}_1.fq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_1_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_1_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1}
        """

rule fastqc_2:
    input:
        r2 = FASTQ_DIR + "/{sample}_2.fq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_2_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_2_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r2}
        """

# Alternative rules for .fastq.gz extension with R1/R2
rule fastqc_r1_alt:
    input:
        r1 = FASTQ_DIR + "/{sample}_R1.fastq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_R1_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_R1_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1}
        """

rule fastqc_r2_alt:
    input:
        r2 = FASTQ_DIR + "/{sample}_R2.fastq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_R2_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_R2_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r2}
        """

# Alternative rules for .fastq.gz extension with _1/_2
rule fastqc_1_alt:
    input:
        r1 = FASTQ_DIR + "/{sample}_1.fastq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_1_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_1_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1}
        """

rule fastqc_2_alt:
    input:
        r2 = FASTQ_DIR + "/{sample}_2.fastq.gz"
    output:
        html = f"{FASTQC_DIR}/{{sample}}_2_fastqc.html",
        zip = f"{FASTQC_DIR}/{{sample}}_2_fastqc.zip"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r2}
        """

# Run MultiQC to aggregate all FastQC results
rule multiqc:
    input:
        # Collect all FastQC outputs for both R1 and R2
        lambda wildcards: expand(f"{FASTQC_DIR}/{{sample}}_{{read}}_fastqc.zip", 
                                sample=SAMPLES, 
                                read=['R1', 'R2'] if any('_R1.' in f for f in r1_files) else ['1', '2'])
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

