##########
#
# Author: F. MASSIP
# date: 12/14/2021
#
# Pipeline to download genomes of different bacteria, compare them, compute the MLD and A for different geographic locations
#

configfile:"config.yml"

SPECIES=config["SPECIES1"]+ config["SPECIES2"]

def get_file_sp(s,mem):
	count_list=[]
	if mem == "l" or mem == "b":
		myfile = config["DATA_DIR_FAST"]+"/processed/"+s+"-country_file_low_mem.txt"
		with open(myfile) as f:
			for i in f.readlines():
				count_list.append(i.rstrip())
	elif mem == "h" or mem == "b" :
		myfileHM = config["DATA_DIR_FAST"]+"/processed/"+s+"-country_file_high_mem.txt"
		with open(myfileHM) as fHM:
			for i in fHM.readlines():
				count_list.append(i.rstrip())
	return count_list


ALL_COUNT1=get_file_sp(config["SPECIES1"],"b")
ALL_COUNT2=get_file_sp(config["SPECIES2"],"b")
COUNTRIES_HM1=get_file_sp(config["SPECIES1"],"h")
COUNTRIES_HM2=get_file_sp(config["SPECIES2"],"h")
COUNTRIES_LM1=get_file_sp(config["SPECIES1"],"l")
COUNTRIES_LM2=get_file_sp(config["SPECIES2"],"l")


rule all:
	input:
#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{COUNTRY1}_{SPECIES2}_{COUNTRY2}.mum",SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES,COUNTRY2=COUNTRIES),
#		expand("{DATA_DIR}/processed/{SPECIES1}/{SPECIES1}_USA.fasta",SPECIES1=config["SPECIES1"],DATA_DIR=config["DATA_DIR"]),
#		expand("{DATA_DIR}/processed/{SPECIES2}/{SPECIES2}_USA.fasta",SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"])
#		expand("{DATA_DIR}/processed/{SPECIES}/{SPECIES}_{COUNTRY}.fasta",SPECIES=SPECIES,DATA_DIR=config["DATA_DIR"],COUNTRY=COUNTRIES),
#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/list_of_files",DATA_DIR=config["DATA_DIR"],SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.mum",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_LM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.h",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_LM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummerHighMem/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.mum",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_HM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummerHighMem/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.h",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_HM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"])



rule cp:
	input:
		fasta1=expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR"]),
	output:
		temp(expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR_FAST"]))
	resources:
		cp_cores=1
	shell:
		"""cp {input.fasta1} {output}
		sleep 5"""


rule mummer:
	input:
		fasta1=expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR_FAST"]),
		fasta2=expand("{DATA_DIR}/processed/{{SPECIES2}}/{{SPECIES2}}_{{COUNTRY2}}.fasta",DATA_DIR=config["DATA_DIR_FAST"])
#		fasta1=get_file_sp1,
#		fasta2=get_file_sp2
	output:
#		listF=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/list_of_files",DATA_DIR=config["DATA_DIR"]),
		mum=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"]),
		h=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.h",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p ~/HGTnew/data/mummer/{wildcards.SPECIES1}_{wildcards.SPECIES2}
		mummer -maxmatch -F -b -n -l 300 {input.fasta1} {input.fasta2} \
		>{output.mum}

		sed '/^>/ d' {output.mum} |\
		sed 's/.* //' |sort -n  >{output.h}"""
	
rule mummerHighMem:
	input:
		fasta1=expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR_FAST"]),
		fasta2=expand("{DATA_DIR}/processed/{{SPECIES2}}/{{SPECIES2}}_{{COUNTRY2}}.fasta",DATA_DIR=config["DATA_DIR_FAST"])
#		fasta1=get_file_sp1,
#		fasta2=get_file_sp2
	output:
#		listF=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/list_of_files",DATA_DIR=config["DATA_DIR"]),
		mum=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"]),
		h=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.h",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p ~/HGTnew/data/mummerHighMem/{wildcards.SPECIES1}_{wildcards.SPECIES2}
		mummer -maxmatch -F -b -n -l 300 {input.fasta1} {input.fasta2} \
		>{output.mum}

		sed '/^>/ d' {output.mum} |\
		sed 's/.* //' |sort -n |uniq -c  >{output.h}"""

#		mummer -maxmatch -b -n -l 300 -F \
#		{config[DATA_DIR]/processed/{wildcards.SPECIES1}/{wildcards.SPECIES1}_{wildcards.COUNTRY1}.fasta \
#		{config[DATA_DIR]/processed/{wildcards.SPECIES2}/{wildcards.SPECIES2}_{wildcards.COUNTRY2}.fasta \
#		>{output}"""

#rule aggregate:
#	input:
#		expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"])
#	output:
#		expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/list_of_files",DATA_DIR=config["DATA_DIR"])
#	shell:
#		"""ls {config["DATA_DIR"]}/processed/mummer/{wildcards.SPECIES1}}_{wildcards.SPECIES2} >{output}"""
	

#		sed '/^>/ d' ~/HGTnew/data/mummer/$sp1\_$sp2/$sp1\_$i\_$sp2\_$j.mum |\
#		sed 's/.* //' |sort -n \
#		>~/HGTnew/data/mummer/$sp1\_$sp2/$sp1\_$i\_$sp2\_$j.mum.h
#		echo $sp1 $i $p2 $j >>$out"""



#		"""{config[CODE_DIR]}/mummer.sh {wildcards.SPECIES1} {wildcards.SPECIES2}\
#		{input.countriesSp1} {input.countriesSp2} {output.listmum}""" 




