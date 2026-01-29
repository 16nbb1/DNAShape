#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=next
#SBATCH --output=next.o
#SBATCH --error=next.e

# trying to run a simple nextflow:
module load nextflow/25.10.2

# Note, within the sample sheet, one vcf corresponds to one individual. See the /home/nboev/projects/def-sushant/nboev/data/CPTAC-3/execute_SampleSheet.sh for the selection I did
mkdir /home/nboev/scratch/nextflow
mkdir -p /home/nboev/scratch/$1/$2/$3
mkdir -p /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3

nextflow run main.nf \
-c hg38.config \
-work-dir /home/nboev/scratch/nextflow \
--samplesheet $4 \
--output_dir /home/nboev/scratch/$1/$2/$3

# Moving the files I will likely need for long term storage (optional)
cp /home/nboev/scratch/$1/$2/$3/merged_5bpDEL.bed /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3
cp /home/nboev/scratch/$1/$2/$3/merged_5bpINS.bed /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3
cp /home/nboev/scratch/$1/$2/$3/merged_5bpSNV.bed /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3

cp /home/nboev/scratch/$1/$2/$3/MERGED_5bpDELSeqcontext.tsv /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3
cp /home/nboev/scratch/$1/$2/$3/MERGED_5bpINSSeqcontext.tsv /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3
cp /home/nboev/scratch/$1/$2/$3/MERGED_5bpSNVSeqcontext.tsv /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/$3


# Running, 2026/01/29
# sbatch tryingnext.sh CPTAC-3 GBM wgs nfSheet_CPTAC3GBM.tsv
# sbatch tryingnext.sh CPTAC-3 PDA wgs nfSheet_CPTAC3PDA.tsv
# sbatch tryingnext.sh CPTAC-3 LUSC wgs nfSheet_CPTAC3LUSC.tsv


# Note: I need to double check what would happen if the vcfs I pass (like for celegans), doesnt have indels?
# In the case for the Zou paper, I should pass in the "constructed" vcfs that contain both?

