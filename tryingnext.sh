#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=next
#SBATCH --output=next.o
#SBATCH --error=next.e

# trying to run a simple nextflow:
module load nextflow/25.10.2

# recall, the samplesheet should already be filtered for patients with primaries + the right cancer type (ie. i should split up my luad/lusc) + pick one vcf per patient?? maybe based on sample quality
# the memory allocation for this job should be super low, but the time needs to span the amount of time we require for all the jobs to finish

nextflow run main.nf \
-c hg38.config \
--samplesheet CPTAC3_BrainWGS.tsv \
--output_dir /home/nboev/scratch/CPTAC-3/BrainCancer/wgs \
-resume

# sbatch tryingnext.sh


# NOTE!! Since I changed the kmer strategy--> need to create a sample-level bed --> intersect with gc content beds, then calculate the kmer counts

# Note: I need to double check what would happen if the vcfs I pass (like for celegans), doesnt have indels?
# In the case for the Zou paper, I should pass in the "constructed" vcfs that contain both?

