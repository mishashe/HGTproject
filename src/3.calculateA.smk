##########
#
# Author: F. MASSIP
# date: 12/14/2021
#

configfile:"config.yml"

rule all:
	input:
		expand("{DATA_DIR}/processed/results/{SPECIES1}_{SPECIES2}.RData",DATA_DIR=config["DATA_DIR"],SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"])


rule CalulateA:
	input:
		mumHM=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}",DATA_DIR=config["DATA_DIR"]),
		mum=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}",DATA_DIR=config["DATA_DIR"])
	output:
		Rdata=expand("{DATA_DIR}/processed/results/{{SPECIES1}}_{{SPECIES2}}.RData",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript {config[CODE_DIR]}/calculate_A_from_mummer_files.R {wildcards.SPECIES1} {wildcards.SPECIES2}"""

