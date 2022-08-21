#!/bin/bash
#SBATCH --job-name=Nicky_cellranger_count_20220520
#SBATCH --output=/project2/tays/Nicky/20220520/sbatch_files/Nicky_cellranger_count_sbatch.out
#SBATCH --error=/project2/tays/Nicky/20220520/sbatch_files/Nicky_cellranger_count_sbatch.err
#SBATCH --partition=bigmem2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem=128G

# Directory variable of working Directory
my_dir=/project2/tays/Nicky/20220520/
mkdir -p $my_dir/Alignment
mkdir -p $my_dir/Alignment/cDNA
cd $my_dir/Alignment/cDNA

# # Load cellranger
# module load cellranger/3.0.0
# Use my own cellranger/6.1.1

# Download human reference genome
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

# Count alignment
/project2/tays/Van/cellranger/cellranger-6.1.1/cellranger count --id=count_10X_pbmc \
--fastqs=$my_dir/FastQ/ \
--sample=cDNA \
--transcriptome=/project2/tays/Van/reference_genomes/10X_human_ref_genome/refdata-gex-GRCh38-2020-A
