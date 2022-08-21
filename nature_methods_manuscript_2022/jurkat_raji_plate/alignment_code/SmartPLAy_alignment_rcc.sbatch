#!/bin/bash
#SBATCH --job-name=Van_SmartPLAy_alignment
#SBATCH --output=/project2/tays/Van/20191001/sbatch_files/Van_SmartPLAy_alignment_sbatch.out
#SBATCH --error=/project2/tays/Van/20191001/sbatch_files/Van_SmartPLAy_alignment_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G

#Directory variable of working Directory
my_dir=/project2/tays/Van/20191001

# Directory for storing summary files
mkdir -p $my_dir/Alignment

#Change directory**
cd $my_dir/Alignment

# Load java**
module load java

java -jar $my_dir/PLA_alignment.jar ReadAlignmentSmartSeq \
R1_LIST=$my_dir/R1List.csv \
O=$my_dir/Alignment/ReadAlignmentSmartSeq_out.txt.gz \
AB_BC_LIST=$my_dir/barcode_cocktail_5.2.csv \
SUMMARY=$my_dir/Alignment/ReadAlignmentSmartSeq_summary.txt \
HEADER=TRUE \
BASE_QUALITY=10 \
NUM_BELOW_BASE_QUALITY=1 \
BASE_QUALITY_START=1

# UMI merging
java -jar $my_dir/PLA_alignment.jar UMIMerging \
I=$my_dir/Alignment/ReadAlignmentSmartSeq_out.txt.gz \
O=$my_dir/Alignment/UMIMerging_out.txt.gz \
SUMMARY=$my_dir/Alignment/UMIMerging_summary.txt

# Digital count matrix
java -jar $my_dir/PLA_alignment.jar DigitalCount \
I=$my_dir/Alignment/UMIMerging_out.txt.gz \
O=$my_dir/Alignment/Jurkat_Raji_plate_count_matrix.txt.gz \
CELL_BC_LIST=NONE \
DUPLICATE_EXPORT=$my_dir/Alignment/DigitalCountduplicate_export.txt.gz \
REMOVE_DUPLICATE=TRUE \
SUMMARY=$my_dir/Alignment/DigitalCount_summary.txt
