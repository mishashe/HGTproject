##########
#
# Author: F. MASSIP
# date: 12/14/2021
#
# Pipeline to download genomes of different bacteria, compare them, compute the MLD and A for different geographic locations
#

configfile:"config.yml"

SPECIES=[config["SPECIES1"],config["SPECIES2"]]

COUNTRIES=["USA","France"]

rule all:
	input:
#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{COUNTRY1}_{SPECIES2}_{COUNTRY2}.mum",SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES,COUNTRY2=COUNTRIES),
#		expand("{DATA_DIR}/processed/{SPECIES1}/{SPECIES1}_USA.fasta",SPECIES1=config["SPECIES1"],DATA_DIR=config["DATA_DIR"]),
#		expand("{DATA_DIR}/processed/{SPECIES2}/{SPECIES2}_USA.fasta",SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"])
#		expand("{DATA_DIR}/processed/{SPECIES}/{SPECIES}_{COUNTRY}.fasta",SPECIES=SPECIES,DATA_DIR=config["DATA_DIR"],COUNTRY=COUNTRIES),
i#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/list_files.txt",SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"]),
		expand("{DATA_DIR}/processed/{SPECIES}_country_file.txt",DATA_DIR=config["DATA_DIR"],SPECIES=config["SPECIES1"]),
		expand("{DATA_DIR}/processed/{SPECIES}_country_file.txt",DATA_DIR=config["DATA_DIR"],SPECIES=config["SPECIES2"])



## Another try to construct a paralellizable rule
def listFile(wildcards):
	checkpoint_output = checkpoints.DownloadSP.get(**wildcards).output[0]
	return expand("{DATA_DIR}/{SPECIES}/{filename}.fasta",DATA_DIR=config["DATA_DIR"],
		SPECIES=wildcards.SPECIES,
		i=glob_wildcards(os.path.join(checkpoint_output, "{filename}.fasta")).{filename})

####
####  Stopped here: I need to find a way to have inputs and outputs that are list of samples from different countries....
####  Maybe I need to output the list of countries in the output file of the Download rule?


rule mummer:
	input:
#		countriesSp1=expand("{DATA_DIR}/processed/{{SPECIES1}}_country_file.txt",DATA_DIR=config["DATA_DIR"]),
#		countriesSp2=expand("{DATA_DIR}/processed/{{SPECIES2}}_country_file.txt",DATA_DIR=config["DATA_DIR"])
		listFile
	output:
		listmum=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/list_files.txt",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p ~/HGTnew/data/mummer/$sp1\_$sp2/
		echo "" >$out
		mummer -maxmatch -b -n -l 300 -F \
		~/HGTnew/data/processed/$sp1/$sp1\_$i\.fasta \
		~/HGTnew/data/processed/$sp2/$sp2\_$j\.fasta \
		>~/HGTnew/data/mummer/$sp1\_$sp2/$sp1\_$i\_$sp2\_$j.mum
		sed '/^>/ d' ~/HGTnew/data/mummer/$sp1\_$sp2/$sp1\_$i\_$sp2\_$j.mum |\
		sed 's/.* //' |sort -n \
		>~/HGTnew/data/mummer/$sp1\_$sp2/$sp1\_$i\_$sp2\_$j.mum.h
		echo $sp1 $i $p2 $j >>$out"""



#		"""{config[CODE_DIR]}/mummer.sh {wildcards.SPECIES1} {wildcards.SPECIES2}\
#		{input.countriesSp1} {input.countriesSp2} {output.listmum}""" 




