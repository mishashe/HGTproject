#!/bin/bash
#SBATCH -J ControlJob
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH --output /cluster/CBIO/data1/data3/fmassip/HGT/ProjectMisha/HGTnew/HGTproject/src/logs/%x-%j.out
#SBATCH --error /cluster/CBIO/data1/data3/fmassip/HGT/ProjectMisha/HGTnew/HGTproject/src/logs/%x-%j.err
eval "$(conda shell.bash hook)"

# activate conda 
conda activate snakemake

# make things fail on errors
set -o nounset
set -o errexit
set -x


export LOGDIR=~/HGTnew/HGTproject/src/logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
### run your commands here!

mkdir -p $LOGDIR


snakemake -s ~/HGTnew/HGTproject/src/1.download.smk -n --unlock

snakemake -s ~/HGTnew/HGTproject/src/1.download.smk \
           --use-conda \
	   --cluster-config config_sge.yml \
	   --cluster "sbatch -N 1 -c 1 -J DL  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
	   --jobs 1 \
	   --rerun-incomplete 

for i in `cat Species_list`
do	
	for j in `cat Species_list`
	do
		if [ $i != $j ]
			then
			my_list=$(echo $i $j | xargs -n1 | sort | xargs)
			cat config_min.yml >config.yml
			echo $my_list |sed -r 's/(.*) (.*)/SPECIES1 : \1\nSPECIES2 : \2/' >>config.yml
			echo $my_list |sed -r 's/(.*) (.*)/SPECIES1 : \1\nSPECIES2 : \2/' >>test

			snakemake -s ~/HGTnew/HGTproject/src/1.8.mash.smk -n --unlock
			snakemake -s ~/HGTnew/HGTproject/src/1.8.mash.smk \
		           --use-conda \
			   --cluster-config config_sge.yml \
			   --cluster "sbatch -N 1 -c 1 -J Mash  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
			   --jobs 10 \
			   --rerun-incomplete \
			   --latency-wait 30


			snakemake -s ~/HGTnew/HGTproject/src/2.compare_genome.smk -n --unlock
				
			snakemake -s ~/HGTnew/HGTproject/src/2.compare_genome.smk \
		           --use-conda \
			   --cluster-config config_sge.yml \
			   --cluster "sbatch -N 1 -c 1 -J Mum  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
			   --jobs 40 \
			   --rerun-incomplete \
			   --latency-wait 60
#			   --resources cp_cores=10 \
			
			snakemake -s ~/HGTnew/HGTproject/src/3.calculateA.smk -n --unlock
			snakemake -s ~/HGTnew/HGTproject/src/3.calculateA.smk \
		           --use-conda \
			   --cluster-config config_sge.yml \
			   --cluster "sbatch -N 1 -c 1 -J CalcA  -o $LOGDIR/%j.log -t {cluster.time} --mem {cluster.mem}" \
			   --jobs 1 \
			   --rerun-incomplete \
			   --latency-wait 30
		fi
	done
done

