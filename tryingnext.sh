#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=00:20:00
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
-c hg38.config \
--input /home/nboev/projects/def-sushant/nboev/preprocess/CPTAC-3/BrainCancer/wgs_alt/aa44f64c-570f-464d-9750-0e1eeede6f45/CPTAC-3.90e5a56c-003c-45b1-9773-e6efc7a51cec.wgs.tumor_normal.gatk4_mutect2.raw_somatic_mutation.vcf \
--output_dir /home/nboev/scratch/CPTAC-3/BrainCancer/wgs/aa44f64c-570f-464d-9750-0e1eeede6f45 \
-resume

# sbatch tryingnext.sh

# RECALL: I WANT TO ADD THE STEP WHERE I OUTPUT THE BED FILES!! VCF TO BED!!! SEE THE

# RECALL: I NEED TO PUT THE PARAMS IN A CONFIG FILE!!!
# RECALL: I NEED TO DOWNLOAD ALL THE PHASTCON DATAINTO THE DB SITE OF SNPEFF
# RECALL: I WILL NEED TO DOWNLOAD SNPEFF DATA FOR HG19 + CELEGANS!!
# RECALL: I DONT NEED TO OUTPUT/SAVE ALL THE OUTPUTS? NEED TO FIGURE OUT WCHICH CAN BE DELETED FROM THE OUTPUTTING DIR (IN SCRATCH?)
# Note: I need to double check what would happen if the vcfs I pass (like for celegans), doesnt have indels?
# In the case for the Zou paper, I should pass in the "constructed" vcfs that contain both?

# to make a config file, i'll need diff ones based on the genome im working with:
# nextflow run main.nf -c base.config, where i can change the config based on the params i need

# possible solution:
# process {publishDir = [path: "results",mode: "copy",cleanup: true]}
# nextflow run main.nf \--input samples/*.vcf \-cleanup \-resume --> essentially add in the cleanup!!but this requires better outputting

# next step is to actually add in the filtering
# for tomorrow, try to create a new manifest file with the actual vcfs without annotations
# in theory, I should: get vcfs > annotate with snpeff (still getting it to work) > dna shape
# i should restructure it so i have original files in the scratch?


# should also add in the part with making a bed file for each sample?--> then intersect with basic files like promoters, enhancers, gene-links? --> maybe keep it sample-level
# keep the merging separate!! doesnt need to be a nextflow

# recall: it needs to also work for the volkova, zou, tcga etc.
# thinking about splitting up the indels/snvs for the processing, maybe for indels, only keeping the means/stds for the wts



