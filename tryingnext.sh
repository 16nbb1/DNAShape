#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=next
#SBATCH --output=next.o
#SBATCH --error=next.e

# trying to run a simple nextflow:
module load nextflow/25.10.2

# had to remove the .gz part, in theory i could add that back in and put it here instead
#names=($(cat ./UterineWGS.txt))
#file=${names[${SLURM_ARRAY_TASK_ID}]}

nextflow run main.nf \
--input /home/nboev/projects/def-sushant/nboev/preprocess/CPTAC-3/BrainCancer/wgs_alt/aa44f64c-570f-464d-9750-0e1eeede6f45/CPTAC-3.90e5a56c-003c-45b1-9773-e6efc7a51cec.wgs.tumor_normal.gatk4_mutect2.raw_somatic_mutation.vcf \
--output_dir /home/nboev/scratch/CPTAC-3/BrainCancer/wgs/aa44f64c-570f-464d-9750-0e1eeede6f45 \
-resume

# I SHOULD SPLIT UP THE SIGPROFILER COMAND AND/OR PASS THE SUBDIR THROUGH!!!
# ex. sigprofiler for indels and one for snvs--> ex. should i put them in one process or should is plit them up?
# alternatively, I could... actually add an argument to the sigprofiler itself to help determine which signatures to use

# RECALL: I NEED TO PUT THE PARAMS IN A CONFIG FILE!!!

# RECALL: I NEED TO DOWNLOAD ALL THE PHASTCON DATAINTO THE DB SITE OF SNPEFF
# good progress, however, I need to delete the extra vcfs that I dont need from the my output dirs (within cancer etc.)
# cant figure out how to delete the intermediate file?
# its not working within the process bc its the input? (ex. rm ${vcf) or other variaitons of this)
# however now i cant figure out how to call it in the next step where i try to run sigprofielr?
# ideally, i could delete teh intermediate file then use this dir for sigprofiler input
# note, i will need to make another dir fo the sigprofiler output


# next step is to actually add in the filtering
# for tomorrow, try to create a new manifest file with the actual vcfs without annotations
# in theory, I should: get vcfs > annotate with snpeff (still getting it to work) > dna shape
# i should restructure it so i have original files in the scratch?


# in the utils, write different functions that are going to be commonly used--> see: mutcontextSNVINDEL.py
#ex: load_data (maybe need to change the name because i have versions with req), save_data, reverse_complement, safe_to_array

# should also add in the part with making a bed file for each sample?--> then intersect with basic files like promoters, enhancers, gene-links? --> maybe keep it sample-level
# keep the merging separate!! doesnt need to be a nextflow

# recall: it needs to also work for the volkova, zou, tcga etc.
# thinking about splitting up the indels/snvs for the processing, maybe for indels, only keeping the means/stds for the wts


# trying to add this to git: MAKE SURE NO VCF INFO IS VISIBLE!! + THE WORK DIR IS CLEANED UP!!
#[nboev@l1.nibi scripts]$ git config --global user.name "Nadejda Boev"
#[nboev@l1.nibi scripts]$ git config --global user.email "nb.boev@mail.utoronto.ca"
#[nboev@l1.nibi scripts]$ git config --global init.defaultBranch main

#git --version
git init
