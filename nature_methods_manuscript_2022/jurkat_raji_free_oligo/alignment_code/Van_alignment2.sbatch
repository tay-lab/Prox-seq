#!/bin/bash
#SBATCH --job-name=Van_SmartPLA_alignment2
#SBATCH --output=/project2/tays/Van/20210722/sbatch_files/Van_SmartPLA_alignment2_sbatch.out
#SBATCH --error=/project2/tays/Van/20210722/sbatch_files/Van_SmartPLA_alignment2_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G

######################
# Use new PLA_alignment.jar with BASE_QUALITY=20 and NUM_BELOW_BASE_QUALITY=1
######################

#Directory variable of working Directory
my_dir=/project2/tays/Van/20210722

# Load java**
module load java

mkdir -p $my_dir/Alignment
mkdir -p $my_dir/Alignment/Van2

# Alignment Van's single cell samples ====================
java -jar $my_dir/PLA_alignment.jar ReadAlignmentSmartSeq \
R1_LIST=$my_dir/R1List_Van.csv \
O=$my_dir/Alignment/Van2/ReadAlignmentSmartSeq_out.txt.gz \
AB_BC_LIST=$my_dir/210518_PLA_cocktail_16_dataanalysis_alignment.csv \
SUMMARY=$my_dir/Alignment/Van2/ReadAlignmentSmartSeq_summary.txt \
HEADER=TRUE \
BASE_QUALITY=20 \
NUM_BELOW_BASE_QUALITY=1 \
BASE_QUALITY_START=1

# UMI merging
java -jar $my_dir/PLA_alignment.jar UMIMerging \
I=$my_dir/Alignment/Van2/ReadAlignmentSmartSeq_out.txt.gz \
O=$my_dir/Alignment/Van2/UMIMerging_out.txt.gz \
SUMMARY=$my_dir/Alignment/Van2/UMIMerging_summary.txt

# Digital count matrix
java -jar $my_dir/PLA_alignment.jar DigitalCount \
I=$my_dir/Alignment/Van2/UMIMerging_out.txt.gz \
O=$my_dir/Alignment/Van2/Van2_count_matrix.txt.gz \
CELL_BC_LIST=NONE \
DUPLICATE_EXPORT=$my_dir/Alignment/Van2/DigitalCountduplicate_export.txt.gz \
REMOVE_DUPLICATE=TRUE \
SUMMARY=$my_dir/Alignment/Van2/DigitalCount_summary.txt
