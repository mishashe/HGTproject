##########
#
# Author: F. MASSIP
# date: 12/14/2021
#
# Pipeline to download genomes of different bacteria, compare them, compute the MLD and A for different geographic locations
#

configfile:"config.yml"


rule all:
	input:
		expand("{DATA_DIR}/processed/mummer/{SPECIES1}_{SPECIES2}.mum",SPECIES1=config["SPECIES1"],SPECIES2=config["SPECIES2"],DATA_DIR=config["DATA_DIR"])



rule Download:
	input:
		csvSP1=expand("{DATA_DIR}/external/{{SPECIES1}}.csv",DATA_DIR=config["DATA_DIR"]),
		csvSP2=expand("{DATA_DIR}/external/{{SPECIES2}}.csv",DATA_DIR=config["DATA_DIR"]),
	output:
		expand("{DATA_DIR}/processed/{{SPECIES1}}/{{SPECIES1}}_USA.fasta",DATA_DIR=config["DATA_DIR"]),
		expand("{DATA_DIR}/processed/{{SPECIES2}}/{{SPECIES2}}_USA.fasta",DATA_DIR=config["DATA_DIR"])

	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript {config[CODE_DIR]}/downloadGenBank.R {wildcards.SPECIES1}
		Rscript {config[CODE_DIR]}/downloadGenBank.R {wildcards.SPECIES2}"""

####
####  Stopped here: I need to find a way to have inputs and outputs that are list of samples from different countries....
####  Maybe I need to output the list of countries in the output file of the Download rule?


rule mummer:
	input:
		data=config["QTL_DIR"] +"data_dir/"+ config["EXP_FILE"]+".tsv",
		cov=config["DATA_DIR"] + "clinical/fixed_clinical.rda",
#		geno=config["DATA_DIR"] + "fgiorgi_project_folder/LungCancerCure/shared/genotypes/fixed_genotypes.rda"
		geno=config["DATA_DIR"] + "phased_and_imputed_v3/combined_batches/final_genotypes.Rda"

	output:
		Rdat=config["QTL_DIR"]  + "results_on_"+config["EXP_FILE"]+"/prepared_for_peer-smoking-grad_and_PacksYear.RData"
	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript  {config[QTL_DIR]}scripts/prepare_peer.R {input.data} {input.geno} {input.cov} {output.Rdat}""" 

rule peer:
	input:
		config["QTL_DIR"]  + "results_on_"+config["EXP_FILE"]+"/prepared_for_peer-smoking-grad_and_PacksYear.RData"
	output:
		pdf=expand("{QTL_DIR}results_on_{EXP_FILE}/{PEER_DIR}{{QQNORM}}/{{NASAL}}/peer_whole_data_cov_{{NB_COV}}.pdf",QTL_DIR=config["QTL_DIR"],PEER_DIR=config["PEER_DIR"],EXP_FILE=config["EXP_FILE"]),
		residuals=expand("{QTL_DIR}results_on_{EXP_FILE}/{PEER_DIR}{{QQNORM}}/{{NASAL}}/residuals_peer_all_samples_cov_{{NB_COV}}.txt",QTL_DIR=config["QTL_DIR"],PEER_DIR=config["PEER_DIR"],EXP_FILE=config["EXP_FILE"])
	params:
		outdir=config["QTL_DIR"]+config["PEER_DIR"]
	conda:
		config["CONDA_FILE"]
	shell:
		"""mkdir -p {params.outdir}/{wildcards.QQNORM};
		Rscript  {config[QTL_DIR]}scripts/peer.R {input} {output.residuals} {output.pdf} {wildcards.NB_COV} {wildcards.QQNORM} {wildcards.NASAL} 2>logs/peer_log_{wildcards.NB_COV}_{wildcards.QQNORM}""" 

#prepare the genotype data in an R object -- Remove colinear snps or snps w/ too low MAF in the pop, and format it for matrixeQTL

rule prep_geno_for_eQTL:
	input:
#		geno=config["INPUT"] + "data_dir/all_genotypes.vcf_reformatted",
#		geno=config["DATA_DIR"] + "fgiorgi_project_folder/LungCancerCure/shared/genotypes/fixed_genotypes.rda",
#		snpMap=config["DATA_DIR"] + "fgiorgi_project_folder/LungCancerCure/lists/snploc_GRCh37.p13.txt"
		geno=config["DATA_DIR"] + "phased_and_imputed_v3/combined_batches/final_genotypes.Rda",
		snpMap=config["DATA_DIR"] + "phased_and_imputed_v3/combined_batches/SnpMap.txt"
#		ID=config["QTL_DIR"] + "ID_equivalent.txt"
	output:
		config["QTL_DIR"] + "preprocessed_genotypes_clinical_maf_sup0.02.RData"
	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript {config[QTL_DIR]}scripts/create_genotype-clinical_object.R  {input.geno} {input.snpMap}"""

# Conduct the eQTL analysis

rule eQTL:
	input:
		geno=config["QTL_DIR"] + "preprocessed_genotypes_clinical_maf_sup0.02.RData",
		pheno=config["QTL_DIR"]  + "results_on_"+config["EXP_FILE"]+ "/prepared_for_peer-smoking-grad_and_PacksYear.RData",
		residuals=expand("{QTL_DIR}results_on_{EXP_FILE}/{PEER_DIR}{{QQNORM}}/{{NASAL}}/residuals_peer_all_samples_cov_{{NB_COV}}.txt",QTL_DIR=config["QTL_DIR"],PEER_DIR=config["PEER_DIR"],EXP_FILE=config["EXP_FILE"])
		
	output:
		pvals=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"]),
		qqplot=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_qqplot.pdf",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"])	
	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript {config[QTL_DIR]}scripts/eQTL_after_peer.R {input.geno} {input.pheno}  {wildcards.NB_COV} {wildcards.COV_TO_TEST} {wildcards.PVAL}  {wildcards.TYPE_PHENO} {wildcards.QQNORM} {config[QTL_DIR]}results_on_{config[EXP_FILE]}/res_using_peer_cov-{wildcards.NB_COV}_pval-{wildcards.PVAL}/{wildcards.TYPE_PHENO}-{wildcards.MODEL}-{wildcards.COV_TO_TEST}-{wildcards.QQNORM}-{wildcards.NASAL} {wildcards.MODEL} {wildcards.NASAL} {input.residuals}"""


rule genotyped_QTLs:
	input:
		pvals=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"]),
		snps=config["DATA_DIR"] + "phased_and_imputed_v2/all_chrs_genotyped_SNPs.txt"
	output:
		pvals_geno=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis_geno",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""python3.6 {config[QTL_DIR]}scripts/select_genotyped_SNPs.py --SNP {input.snps} -e {input.pvals} >{output.pvals_geno}"""

rule GeneName:
	input:
		QTLs=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"]),
		QTLs_geno=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis_geno",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"])
	output:
		annot=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis_w_gene_names",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"]),
		annot_geno=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis_geno_w_gene_names",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"])
	conda:
		config["CONDA_FILE"]
	shell:
		"""{config[QTL_DIR]}scripts/find_gene_name.sh {input.QTLs} {wildcards.COV_TO_TEST}-{wildcards.PVAL}-{wildcards.QQNORM}-{wildcards.NB_COV}-{wildcards.NASAL}
		{config[QTL_DIR]}scripts/find_gene_name.sh {input.QTLs_geno} {wildcards.COV_TO_TEST}-{wildcards.PVAL}-{wildcards.QQNORM}-{wildcards.NB_COV}-{wildcards.NASAL}"""

rule PostQTL:
	input:
		geno=config["QTL_DIR"] + "preprocessed_genotypes_clinical_maf_sup0.02.RData",
		pheno=config["QTL_DIR"]  +"results_on_"+config["EXP_FILE"]+ "/prepared_for_peer-smoking-grad_and_PacksYear.RData",
		residuals=expand("{QTL_DIR}results_on_{EXP_FILE}/{PEER_DIR}{{QQNORM}}/{{NASAL}}/residuals_peer_all_samples_cov_{{NB_COV}}.txt",QTL_DIR=config["QTL_DIR"],PEER_DIR=config["PEER_DIR"],EXP_FILE=config["EXP_FILE"]),
		QTLs=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis_w_gene_names",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"])
	output:
		pdf=expand("{QTL_DIR}results_on_{EXP_FILE}/res_using_peer_cov-{{NB_COV}}_pval-{{PVAL}}/{{TYPE_PHENO}}-{{MODEL}}-{{COV_TO_TEST}}-{{QQNORM}}-{{NASAL}}/Matrix_eQTL_analysis.txt_cis_w_gene_names_all_samples_boxplot_Cancer.pdf",QTL_DIR=config["QTL_DIR"],EXP_FILE=config["EXP_FILE"]),
	conda:
		config["CONDA_FILE"]
	shell:
		"""Rscript {config[QTL_DIR]}scripts/analysis_post_QTL.R {input.geno} {input.pheno}  {wildcards.NB_COV} {wildcards.COV_TO_TEST} {wildcards.PVAL}  {wildcards.TYPE_PHENO} {wildcards.QQNORM} {input.QTLs} {wildcards.MODEL} {wildcards.NASAL} {input.residuals}"""



