#!/bin/bash
#SBATCH --job-name=Van_cellranger_count_downsample_20210831
#SBATCH --output=/project2/tays/Van/20210831/sbatch_files/cellranger_count_downsample_sbatch.out
#SBATCH --error=/project2/tays/Van/20210831/sbatch_files/cellranger_count_downsample_sbatch.err
#SBATCH --partition=bigmem2
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem=128G

# Directory variable of working Directory
my_dir=/project2/tays/Van/20210831
mkdir -p "$my_dir"/downsample
mkdir -p "$my_dir"/downsample/cDNA
mkdir -p "$my_dir"/downsample/cDNA/FastQ
cd "$my_dir"/downsample/cDNA

# # Load cellranger
# module load cellranger/3.0.0
# Use my own cellranger/6.1.1
cellranger=/project2/tays/Van/cellranger/cellranger-6.1.1/cellranger

# Download human reference genome
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

# Download seqkit for read downsampling
# wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
seqkit=/project2/tays/Van/seqkit/seqkit
rm -f "$my_dir"/downsample/cDNA/seqkit.out
rm -f "$my_dir"/downsample/cDNA/seqkit.err

for i in 0.1 0.2 0.4 0.6 0.8; do
  # Reformat i: replace . with _ (cellranger doesn't handle . well)
  j="${i/./_}"

  # Downsample R1
  # echo $(date) >> "$my_dir"/downsample/cDNA/seqkit.out
  echo $(date) >> "$my_dir"/downsample/cDNA/seqkit.err
  $seqkit sample "$my_dir"/FastQ/cDNA_S1_R1_001.fastq.gz \
    -p "$i" -s 0 --two-pass \
    -o "$my_dir"/downsample/cDNA/FastQ/cDNA_"$j"_S1_R1_001.fastq.gz \
    >> "$my_dir"/downsample/cDNA/seqkit.out 2>> "$my_dir"/downsample/cDNA/seqkit.err

  # Downsample R2
  # echo $(date) >> "$my_dir"/downsample/cDNA/seqkit.out
  echo $(date) >> "$my_dir"/downsample/cDNA/seqkit.err
  $seqkit sample "$my_dir"/FastQ/cDNA_S1_R2_001.fastq.gz \
    -p "$i" -s 0 --two-pass \
    -o "$my_dir"/downsample/cDNA/FastQ/cDNA_"$j"_S1_R2_001.fastq.gz \
    >> "$my_dir"/downsample/cDNA/seqkit.out 2>> "$my_dir"/downsample/cDNA/seqkit.err
  # zcat "$my_dir"/FastQ/cDNA_S1_R2_001.fastq.gz \
  #   | $seqkit sample -p "$i" -s 0 -o "$my_dir"/downsample/cDNA/FastQ/cDNA_"$i"_S1_R2_001.fastq.gz \
  #   >> "$my_dir"/downsample/cDNA/seqkit.out 2>> "$my_dir"/downsample/cDNA/seqkit.err

  # Alignment
  $cellranger count \
    --id=count_downsample_"$j" \
    --fastqs="$my_dir"/downsample/cDNA/FastQ \
    --transcriptome=/project2/tays/Van/reference_genomes/10X_human_ref_genome/refdata-gex-GRCh38-2020-A \
    --sample=cDNA_"$j" \
    --expect-cells=10000 \
    --nosecondary --no-bam

  # Remove downsample fastq files
  rm "$my_dir"/downsample/cDNA/FastQ/cDNA_"$j"_S1_R1_001.fastq.gz
  rm "$my_dir"/downsample/cDNA/FastQ/cDNA_"$j"_S1_R2_001.fastq.gz
done
