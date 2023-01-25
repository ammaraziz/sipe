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

Activate conda env
```
conda activate sipe
```

sample
```
snakemake -j {threads} -k
```

### Install

Install conda then mamba. Then run

```
mamba create -n sipe -c bioconda -c conda-forge -k ivar bowtie2 fastp bamtools snakemake
```

### Depends

Install in order:
1. conda/mamba
2. ivar
3. bowtie2
4. samtools
5. bamtools
6. fastp
