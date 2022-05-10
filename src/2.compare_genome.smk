##########
#
# Author: F. MASSIP
# date: 12/14/2021
#

configfile:"config.yml"

SPECIES=config["SPECIES1"]+ config["SPECIES2"]

def get_file_sp(s,mem):
	count_list=[]
	if mem == "l" or mem == "b":
		myfile = config["DATA_DIR"]+"/processed/"+s+"-country_file_low_mem.txt"
		with open(myfile) as f:
			for i in f.readlines():
				count_list.append(i.rstrip())
	if mem == "h" or mem == "b" :
		myfileHM = config["DATA_DIR"]+"/processed/"+s+"-country_file_high_mem.txt"
		with open(myfileHM) as fHM:
			for i in fHM.readlines():
				count_list.append(i.rstrip())
	if mem == "t" :
		myfileHM = config["DATA_DIR"]+"/processed/"+s+"-country_file_test.txt"
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
#		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/list_of_files",DATA_DIR=config["DATA_DIR"],SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.mum",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_LM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.h",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_LM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummerHighMem/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.mum",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_HM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummerHighMem/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.h",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_HM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.mum.h-filtered-mash",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_LM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/mummerHighMem/{SPECIES1}_{SPECIES2}/{SPECIES1}-{COUNTRY1}_{SPECIES2}-{COUNTRY2}.mum.h-filtered-mash",DATA_DIR=config["DATA_DIR"],COUNTRY1=COUNTRIES_HM1,COUNTRY2=ALL_COUNT2,SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"]),
		expand("{DATA_DIR}/processed/toFilterMash/{SPECIES1}/{SPECIES1}-{COUNTRY1}/toFilter.txt",DATA_DIR=config["DATA_DIR"],COUNTRY1=ALL_COUNT1,SPECIES1=config["SPECIES1"]),
		expand("{DATA_DIR}/processed/toFilterMash/{SPECIES2}/{SPECIES2}-{COUNTRY2}/toFilter.txt",DATA_DIR=config["DATA_DIR"],COUNTRY2=ALL_COUNT2,SPECIES2=config["SPECIES2"])


rule mummer:
	input:
		fasta1=expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR"]),
		fasta2=expand("{DATA_DIR}/processed/{{SPECIES2}}/{{SPECIES2}}_{{COUNTRY2}}.fasta",DATA_DIR=config["DATA_DIR"])
	output:
		mum=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"]),
		h=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.h",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p ~/HGTnew/data/mummer/{wildcards.SPECIES1}_{wildcards.SPECIES2}
		mummer -maxmatch -F -b -n -l 100 {input.fasta1} {input.fasta2} \
		>{output.mum}

		sed '/^>/ d' {output.mum} |\
		sed 's/.* //' |sort -n |uniq -c >{output.h}"""
	
rule mummerHighMem:
	input:
		fasta1=expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR"]),
		fasta2=expand("{DATA_DIR}/processed/{{SPECIES2}}/{{SPECIES2}}_{{COUNTRY2}}.fasta",DATA_DIR=config["DATA_DIR"])
	output:
		mum=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"]),
		h=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.h",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p ~/HGTnew/data/mummerHighMem/{wildcards.SPECIES1}_{wildcards.SPECIES2}
		mummer -maxmatch -F -b -n -l 100 {input.fasta1} {input.fasta2} \
		>{output.mum}

		sed '/^>/ d' {output.mum} |\
		sed 's/.* //' |sort -n |uniq -c  >{output.h}"""

rule filterHighMem:
	input:
		mum=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"]),
		filtFile=expand("{DATA_DIR}/processed/toFilterMash/{{SPECIES1}}/{{SPECIES1}}-{{COUNTRY1}}/toFilter.txt",DATA_DIR=config["DATA_DIR"])
	output:
		hfilt=expand("{DATA_DIR}/processed/mummerHighMem/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum.h-filtered-mash",DATA_DIR=config["DATA_DIR"]),

	conda:
		config["CONDA_FILE"]
	shell:
		"""python3 {config[CODE_DIR]}filter_mummer_file_mash.py --species1 {wildcards.SPECIES1}\
		--species2 {wildcards.SPECIES2} \
		--country1 {wildcards.COUNTRY1} --country2 {wildcards.COUNTRY2} \
		-f {input.mum}"""

rule filterlowMem:
	input:
		mum=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum",DATA_DIR=config["DATA_DIR"]),
		filtFile=expand("{DATA_DIR}/processed/toFilterMash/{{SPECIES1}}/{{SPECIES1}}-{{COUNTRY1}}/toFilter.txt",DATA_DIR=config["DATA_DIR"])

	output:
		hfilt=expand("{DATA_DIR}/processed/mummer/{{SPECIES1}}_{{SPECIES2}}/{{SPECIES1}}-{{COUNTRY1}}_{{SPECIES2}}-{{COUNTRY2}}.mum.h-filtered-mash",DATA_DIR=config["DATA_DIR"]),
	conda:
		config["CONDA_FILE"]
	shell:
		"""python3 {config[CODE_DIR]}filter_mummer_file_mash.py --species1 {wildcards.SPECIES1} \
		--species2 {wildcards.SPECIES2} \
		--country1 {wildcards.COUNTRY1} --country2 {wildcards.COUNTRY2} \
		-f {input.mum}"""


