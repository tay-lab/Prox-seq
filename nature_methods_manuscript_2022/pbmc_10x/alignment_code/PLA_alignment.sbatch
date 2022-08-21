#!/bin/bash
#SBATCH --job-name=10x_PLA_alignment
#SBATCH --output=/project2/tays/Van/20210831/sbatch_files/10x_PLA_alignment_sbatch.out
#SBATCH --error=/project2/tays/Van/20210831/sbatch_files/10x_PLA_alignment_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G

#Directory variable of working Directory
my_dir=/project2/tays/Van/20210831

# Load java**
module load java

mkdir -p $my_dir/Alignment
mkdir -p $my_dir/Alignment/PLA

# Alignment
java -jar $my_dir/PLA_alignment.jar ReadAlignment10x \
R1=$my_dir/FastQ_PLA_combined/PLA_S1_R1_001.fastq.gz \
R2=$my_dir/FastQ_PLA_combined/PLA_S1_R2_001.fastq.gz \
O=$my_dir/Alignment/PLA/ReadAlignment10x_out.txt.gz \
AB_BC_LIST=$my_dir/cocktail_2.6_alignment.csv \
SUMMARY=$my_dir/Alignment/PLA/ReadAlignment10x_summary.txt \
HEADER=TRUE \
BASE_QUALITY=10 \
NUM_BELOW_BASE_QUALITY=1 \
BASE_QUALITY_START=1

# Cell barcode correction
java -jar $my_dir/PLA_alignment.jar CellBarcodeCorrection \
I=$my_dir/Alignment/PLA/ReadAlignment10x_out.txt.gz \
O=$my_dir/Alignment/PLA/CellBarcodeCorrection_out.txt.gz \
MODE=10X \
CELL_BC_LIST=$my_dir/Alignment/cDNA/count_10X_pbmc/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
SUMMARY=$my_dir/Alignment/PLA/CellBarcodeCorrection_summary.txt \
SUFFIX=-1 \
HEADER=FALSE

# UMI merging
java -jar $my_dir/PLA_alignment.jar UMIMerging \
I=$my_dir/Alignment/PLA/CellBarcodeCorrection_out.txt.gz \
O=$my_dir/Alignment/PLA/UMIMerging_out.txt.gz \
SUMMARY=$my_dir/Alignment/PLA/UMIMerging_summary.txt

# Digital count matrix
java -jar $my_dir/PLA_alignment.jar DigitalCount \
I=$my_dir/Alignment/PLA/UMIMerging_out.txt.gz \
O=$my_dir/Alignment/PLA/10x_PLA_count_matrix.txt.gz \
CELL_BC_LIST=NONE \
DUPLICATE_EXPORT=$my_dir/Alignment/PLA/DigitalCountduplicate_export.txt.gz \
REMOVE_DUPLICATE=TRUE \
SUMMARY=$my_dir/Alignment/PLA/DigitalCount_summary.txt
