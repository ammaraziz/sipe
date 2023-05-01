# Snakemake Illumina Pipeline for Entero - sipe

### Overview

1. fastp 	- trim and qc 
2. blast	- blast reads against vp1 reference
3. Bowtie2 	- Align to full genome nr entero reference set 
4. ivar		- create consensus

### Inputs

- config.yaml
	- out_dir: path to output directory
	- reference_wgf: path to reference for alignment
	- reference_blast: path to blast database
	- bed_file: path to primer bed files - not used
	
- samples.tsv
	- columns: forward reverse sample sample_library

### Outputs

1. See `analysis/` folder
	- `qc/{sample}.{pair}.fastq.gz` - qc'ed fastq
	- `qc/{sample}.html`
	- `blast/{sample}.xml` - blast results
	- `align/{sample}/{sample}.bam`
	- `consensus/{sample}/{sample}.consensus.fasta`
	- `depth/{sample}.depth.txt`
	- `depth/{sample}.pdf`

### Install

Install conda/mamba and snakeamke

```
mamba create -c bioconda -c conda-forge snakemake
```

### Run:

Navigate to `sipe/` then run:

```
snakemake -j {cpu} --use-conda
```

Snakemake will handle the dependency installation.

### Depends

Install in order:
1. conda/mamba
2. ivar
3. bowtie2
4. samtools
5. bamtools
6. fastp
