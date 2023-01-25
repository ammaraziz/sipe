# Snakemake Illumina Pipeline for Entero - sipe

### Overview

1. fastp 	- trim and qc 
2. Bowtie2 	- Align to VP1 nr entero reference set using 
3. ivar 	- trim primers
4. ivar		- create consensus
5. mosdepth - generate depth coverage report

### Inputs

- config.yaml
	- samples_path: path to sample.tsv
	- out_dir: path to output directory
	- reference[vp1, wgs]: specify reference set, vp or whole genome
	- bed_file: path to primer bed files
	- library_delimiter['_']: file name delimiter 
	
- samples.tsv
	- columns: forward reverse sample sample_library

### Outputs

TBA

1. See `analysis/` folder
	- `alignments/{sample}/{sample}.bam`
	- `consensus/{sample}/{sample}.consensus.fasta`
	- `depth/{sample}/{sample}.depth.txt`
	- `depth/{sample}/depth.pdf`
	- `qc/{sample}.html`

### Run

```
snakemake -j {threads} -k
```

### Install

Install conda then mamba. Then run

```
tba
```

### Depends

1. conda/mamba
2. bowtie2
3. samtools
4. bamtools
4. fastp
5. mosdepth
7. ivar - install  this first
