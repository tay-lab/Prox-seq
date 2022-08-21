#!/bin/bash
#SBATCH --job-name=Van_PLA_alignment_downsample_20210831
#SBATCH --output=/project2/tays/Van/20210831/sbatch_files/PLA_alignment_downsample_sbatch.out
#SBATCH --error=/project2/tays/Van/20210831/sbatch_files/PLA_alignment_downsample_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G

# Directory variable of working Directory
my_dir=/project2/tays/Van/20210831
mkdir -p "$my_dir"/downsample
mkdir -p "$my_dir"/downsample/PLA
mkdir -p "$my_dir"/downsample/PLA/FastQ
cd "$my_dir"/downsample/PLA

# Load java**
module load java

# Download seqkit for read downsampling
# wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
seqkit=/project2/tays/Van/seqkit/seqkit
rm -f "$my_dir"/downsample/PLA/seqkit.out
rm -f "$my_dir"/downsample/PLA/seqkit.err

for i in 0.1 0.2 0.4 0.6 0.8; do
  # Reformat i: replace . with _ (cellranger doesn't handle . well)
  j="${i/./_}"

  # Make a new folder
  mkdir -p "$my_dir"/downsample/PLA/count_downsample_"$j"

  # Downsample R1
  # echo $(date) >> "$my_dir"/downsample/PLA/seqkit.out
  echo $(date) >> "$my_dir"/downsample/PLA/seqkit.err
  $seqkit sample "$my_dir"/FastQ_PLA_combined/PLA_S1_R1_001.fastq.gz \
    -p "$i" -s 0 --two-pass \
    -o "$my_dir"/downsample/PLA/FastQ/PLA_"$j"_S1_R1_001.fastq.gz \
    >> "$my_dir"/downsample/PLA/seqkit.out 2>> "$my_dir"/downsample/PLA/seqkit.err

  # Downsample R2
  # echo $(date) >> "$my_dir"/downsample/PLA/seqkit.out
  echo $(date) >> "$my_dir"/downsample/PLA/seqkit.err
  $seqkit sample "$my_dir"/FastQ_PLA_combined/PLA_S1_R2_001.fastq.gz \
    -p "$i" -s 0 --two-pass \
    -o "$my_dir"/downsample/PLA/FastQ/PLA_"$j"_S1_R2_001.fastq.gz \
    >> "$my_dir"/downsample/PLA/seqkit.out 2>> "$my_dir"/downsample/PLA/seqkit.err
  # zcat "$my_dir"/FastQ/PLA_S1_R2_001.fastq.gz \
  #   | $seqkit sample -p "$i" -s 0 -o "$my_dir"/downsample/PLA/FastQ/PLA_"$i"_S1_R2_001.fastq.gz \
  #   >> "$my_dir"/downsample/PLA/seqkit.out 2>> "$my_dir"/downsample/PLA/seqkit.err

  # Alignment
  # Alignment
  java -jar "$my_dir"/PLA_alignment.jar ReadAlignment10x \
  R1="$my_dir"/downsample/PLA/FastQ/PLA_"$j"_S1_R1_001.fastq.gz \
  R2="$my_dir"/downsample/PLA/FastQ/PLA_"$j"_S1_R2_001.fastq.gz \
  O="$my_dir"/downsample/PLA/count_downsample_"$j"/ReadAlignment10x_out.txt.gz \
  AB_BC_LIST="$my_dir"/cocktail_2.6_alignment.csv \
  SUMMARY="$my_dir"/downsample/PLA/count_downsample_"$j"/ReadAlignment10x_summary.txt \
  HEADER=TRUE \
  BASE_QUALITY=10 \
  NUM_BELOW_BASE_QUALITY=1 \
  BASE_QUALITY_START=1

  # Cell barcode correction
  java -jar "$my_dir"/PLA_alignment.jar CellBarcodeCorrection \
  I="$my_dir"/downsample/PLA/count_downsample_"$j"/ReadAlignment10x_out.txt.gz \
  O="$my_dir"/downsample/PLA/count_downsample_"$j"/CellBarcodeCorrection_out.txt.gz \
  MODE=10X \
  CELL_BC_LIST="$my_dir"/downsample/cDNA/count_downsample_"$j"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
  SUMMARY="$my_dir"/downsample/PLA/count_downsample_"$j"/CellBarcodeCorrection_summary.txt \
  SUFFIX=-1 \
  HEADER=FALSE

  # UMI merging
  java -jar "$my_dir"/PLA_alignment.jar UMIMerging \
  I="$my_dir"/downsample/PLA/count_downsample_"$j"/CellBarcodeCorrection_out.txt.gz \
  O="$my_dir"/downsample/PLA/count_downsample_"$j"/UMIMerging_out.txt.gz \
  SUMMARY="$my_dir"/downsample/PLA/count_downsample_"$j"/UMIMerging_summary.txt

  # Digital count matrix
  java -jar "$my_dir"/PLA_alignment.jar DigitalCount \
  I="$my_dir"/downsample/PLA/count_downsample_"$j"/UMIMerging_out.txt.gz \
  O="$my_dir"/downsample/PLA/count_downsample_"$j"/10x_PLA_count_matrix.txt.gz \
  CELL_BC_LIST=NONE \
  DUPLICATE_EXPORT="$my_dir"/downsample/PLA/count_downsample_"$j"/DigitalCountduplicate_export.txt.gz \
  REMOVE_DUPLICATE=TRUE \
  SUMMARY="$my_dir"/downsample/PLA/count_downsample_"$j"/DigitalCount_summary.txt

  # Remove downsample fastq files
  rm "$my_dir"/downsample/PLA/FastQ/PLA_"$j"_S1_R1_001.fastq.gz
  rm "$my_dir"/downsample/PLA/FastQ/PLA_"$j"_S1_R2_001.fastq.gz
done
