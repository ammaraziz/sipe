import os
import re
import time
import pandas as pd
from pathlib import Path
from shutil import copyfile
from datetime import datetime

configfile: "config/config.json"
reference = Path(config["reference"])
bed_file = Path(config["bed_file"])
out_dir = Path(config["out_dir"])
input_dir = Path(config["input_dir"])
#samples_path = Path(config["samples_path"])

# get sample names 
#samples = pd.read_table(samples_path, sep="\t")
SAMPLE_NAME, SAMPLE_NUMBER, PAIR = glob_wildcards(input_dir / "{sample_name}_{sample_number}_L001_{pair}_001.fastq.gz")
SAMPLES = list(set([i + "_" + x for i, x in zip(SAMPLE_NAME, SAMPLE_NUMBER)]))
print(input_dir)
rule all:
    input:
        expand(out_dir / "qc/{sample}.html", sample = SAMPLES),
        expand(out_dir / "align/{sample}/{sample}.bam", sample = SAMPLES),
        expand(out_dir / "status/{sample}.split.bam", sample = SAMPLES),
        expand(out_dir / "depth/raw/{sample}.txt", sample = SAMPLES),
        #expand(out_dir / "depth/raw/{sample}.png", sample = SAMPLES),

rule qc:
	message: "Running fastp on {wildcards.sample}"
	input:
		forward = expand(input_dir / "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
		rev = expand(input_dir / "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
	output:
		report = out_dir / "qc/{sample}.html",
		stats = out_dir / "qc/{sample}.stats",
	shell:"""
	fastp -i {input.forward} -I {input.rev} -h {output.report} 2> {output.stats}
	"""

# rule build_index:
# 	message: "Building bowtie2 index"
# 	input:
# 		reference = reference
# 	output:
# 		status = out_dir / "status/bowtie2.index"
# 	params:
# 		prefix = reference.parent / (reference.name + ".index")
# 	shell:"""
# 	bowtie2-build -q {input.reference} {params.prefix}
# 	touch {output}
#	"""

rule align:
	message: "Aligning with bowtie2"
	input:
		forward = expand(input_dir / "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R1"]),
		rev = expand(input_dir / "{{sample}}_L001_{pair}_001.fastq.gz", pair = ["R2"]),
	output:
		bam = out_dir / "align/{sample}/{sample}.bam",
		stats = out_dir / "align/{sample}/{sample}.stats",
	params:
		index = reference,
	threads: 20
	shell:"""
	bowtie2 -q \
	--no-unal \
	-p {threads} \
	-x {params.index} \
	-1 {input.forward} \
	-2 {input.rev} \
	--no-mixed \
	--no-discordant \
	2> {output.stats} \
	| samtools view -bS - \
	| samtools sort -@ {threads} - 1> {output.bam} 2> /dev/null
	samtools index {output.bam}
	"""

rule split_bam:
	message: "Splitting BAM files by chromosome"
	input:
		bam = rules.align.output.bam
	output:
		status = out_dir / "status/{sample}.split.bam"
	params:
		split_prefix = out_dir / "align/{sample}/"
	shell:"""
	bamtools split -in {input.bam} -refPrefix {params.split_prefix} -reference
	touch {output.status}
	"""

rule depth:
	message: "Calculating depth per alignment"
	input:
		bam = rules.align.output.bam
	output:
		depth = out_dir / "depth/raw/{sample}.txt"
	shell:"""
	samtools depth -a {input.bam} -o {output.depth}
	"""

# rule plot_depth:
# 	message: "plotting depth"
# 	input:
# 		depth = rules.depth.output.plot
# 	output:
# 		png = outdir / "depth/raw/{sample}.png"
# 	run:
# 		import .scripts










