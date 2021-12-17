##########
#
# Author: F. MASSIP
# date: 12/14/2021
#
# Pipeline to download genomes of different bacteria, compare them, compute the MLD and A for different geographic locations
#

configfile:"config.yml"

SPECIES=[config["SPECIES1"],config["SPECIES2"]]


rule all:
	input:
#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{COUNTRY1}_{SPECIES2}_{COUNTRY2}.mum",SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES,COUNTRY2=COUNTRIES),
#		expand("{DATA_DIR}/processed/{SPECIES1}/{SPECIES1}_USA.fasta",SPECIES1=config["SPECIES1"],DATA_DIR=config["DATA_DIR"]),
#		expand("{DATA_DIR}/processed/{SPECIES2}/{SPECIES2}_USA.fasta",SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"])
#		expand("{DATA_DIR}/processed/{SPECIES}/{SPECIES}_{COUNTRY}.fasta",SPECIES=SPECIES,DATA_DIR=config["DATA_DIR"],COUNTRY=COUNTRIES),
i#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/list_files.txt",SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"]),
		expand("{DATA_DIR}/processed/{SPECIES}_country_file.txt",DATA_DIR=config["DATA_DIR_FAST"],SPECIES=config["SPECIES1"]),
		expand("{DATA_DIR}/processed/{SPECIES}_country_file.txt",DATA_DIR=config["DATA_DIR_FAST"],SPECIES=config["SPECIES2"])


rule DownloadSP:
	input:
		csvSP1=expand("{DATA_DIR}/external/{{SPECIES}}.csv",DATA_DIR=config["DATA_DIR_FAST"])
	output:
#		fastaSp1=expand("{DATA_DIR}/processed/{{SPECIES}}/{{SPECIES}}_{{COUNTRY}}.fasta",DATA_DIR=config["DATA_DIR"]),
#		countries=expand("{DATA_DIR}/processed/{{SPECIES}}_country_file.txt",DATA_DIR=config["DATA_DIR"]),
		fasta=directory(config["DATA_DIR_FAST"]+"/processed/{SPECIES}")
	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript {config[CODE_DIR]}/downloadGenBank.R {wildcards.SPECIES} {output.countries}"""



