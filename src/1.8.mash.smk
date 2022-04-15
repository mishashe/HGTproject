##########
#
# Author: F. MASSIP
# date: 12/14/2021
#

configfile:"config.yml"

SPECIES=config["SPECIES1"]+ config["SPECIES2"]

def get_countries(s,mem):
	count_list=[]
	if mem == "t" :
		myfileHM = config["DATA_DIR"]+"/processed/"+s+"-country_file_test.txt"
		with open(myfileHM) as fHM:
			for i in fHM.readlines():
				count_list.append(i.rstrip())
	return count_list

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
	if mem == "t1" :
		myfileHM = config["DATA_DIR"]+"/processed/"+s+"-country_file_test.txt"
		with open(myfileHM) as fHM:
			for i in fHM.readlines():
				count_list.append(i.rstrip())

	count_pairs=[]
	for i in range(0,(len(count_list)-1)):
		for j in range(i+1,(len(count_list))):
			count_pairs.append('-'.join(count_list[u] for u in [i,j]))
	return count_pairs

COUNTRY_PAIRS_SP1 = get_file_sp(config["SPECIES1"],"b")
COUNTRY_PAIRS_SP2 = get_file_sp(config["SPECIES2"],"b")

COUNTRIES = get_countries(config["SPECIES1"],"b")

rule all:
	input:
#		expand("{DATA_DIR}/processed/nucmer/{SPECIES1}/{SPECIES1}_{COUNTRIES}.delta",DATA_DIR=config["DATA_DIR"],COUNTRIES=COUNTRY_PAIRS_SP1,SPECIES1=config["SPECIES1"]),
#		expand("{DATA_DIR}/processed/nucmer/{SPECIES1}/{SPECIES1}_{COUNTRIES}.coo",DATA_DIR=config["DATA_DIR"],COUNTRIES=COUNTRY_PAIRS_SP1,SPECIES1=config["SPECIES1"]),
#		expand("{DATA_DIR}/processed/nucmer/{SPECIES2}/{SPECIES2}_{COUNTRIES}.delta",DATA_DIR=config["DATA_DIR"],COUNTRIES=COUNTRY_PAIRS_SP2,SPECIES2=config["SPECIES2"]),
#		expand("{DATA_DIR}/processed/nucmer/{SPECIES2}/{SPECIES2}_{COUNTRIES}.coo",DATA_DIR=config["DATA_DIR"],COUNTRIES=COUNTRY_PAIRS_SP2,SPECIES2=config["SPECIES2"])
		expand("{DATA_DIR}/dist/{SPECIES1}-{COUNTRIES}/dist-all.txt",DATA_DIR=config["DATA_DIR"],COUNTRIES=COUNTRY_PAIRS_SP1,SPECIES1=config["SPECIES1"]),
		expand("{DATA_DIR}/dist/{SPECIES2}-{COUNTRIES}/dist-all.txt",DATA_DIR=config["DATA_DIR"],COUNTRIES=COUNTRY_PAIRS_SP2,SPECIES2=config["SPECIES2"])



rule sketch:
	input:
		fasta1=expand("{DATA_DIR}/processed/{{SPECIES}}/{{SPECIES}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR"]),
	output:
		sketchdir=expand("{DATA_DIR}/tmp/{{SPECIES}}-{{COUNTRY1}}-sketch/",DATA_DIR=config["DATA_DIR"]),

	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p {output.sketchdir}
		demulty_fasta -i {input.fasta1} -o {output.sketchdir}"""

rule mash:
	input:
		sketch1=expand("{DATA_DIR}/tmp/{{SPECIES}}-{{COUNTRY1}}-sketch/",DATA_DIR=config["DATA_DIR"]),
		fasta2=expand("{DATA_DIR}/processed/{{SPECIES}}/{{SPECIES}}_{{COUNTRY1}}.fasta",DATA_DIR=config["DATA_DIR"])
	output:
		dist1=expand("{DATA_DIR}/dist/{{SPECIES}}-{{COUNTRY1}}-{{COUNTRY2}}/dist-all.txt",DATA_DIR=config["DATA_DIR"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p {config[DATA_DIR]}/dist/{wildcards.SPECIES}-{wildcards.COUNTRY1}-{wildcards.COUNTRY2}
		for i in `ls {input.sketch1}/*msh`; \
		do; \
		mash dist $i {input.fasta2} >>{output.dist1} ;\
		done"""

	
