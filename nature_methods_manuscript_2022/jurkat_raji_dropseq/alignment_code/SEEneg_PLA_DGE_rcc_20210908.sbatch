#!/bin/bash
#SBATCH --job-name=Van_SEEneg_PLA_DGE_20210908
#SBATCH --output=/project2/tays/Van/20200104/sbatch_files/Van_SEEneg_PLA_DGE_20210908_sbatch.out
#SBATCH --error=/project2/tays/Van/20200104/sbatch_files/Van_SEEneg_PLA_DGE_20210908_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G

# Directory variable of working Directory
my_dir=/project2/tays/Van/20200104
cd $my_dir/Alignments/

# Make a folder for each sample in the Alignments directory**
mkdir -p SEEneg-PLA

# Load java**
module load java

# # Align raw reads
# java -jar $my_dir/PLA_alignment.jar ReadAlignmentDropSeq \
# R1=$my_dir/FastQ/N703-SEEneg-PLA_S3_R1_001.fastq.gz \
# R2=$my_dir/FastQ/N703-SEEneg-PLA_S3_R2_001.fastq.gz \
# AB_BC_LIST=$my_dir/Alignments/barcode_cocktail_5.3.csv \
# O=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned.txt.gz \
# SUMMARY=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_summary.txt \
# HEADER=TRUE \
# BASE_QUALITY=10 \
# NUM_BELOW_BASE_QUALITY=1 \
# BASE_QUALITY_START=1
#
# # Cell barcode correction
# java -jar $my_dir/PLA_alignment.jar CellBarcodeCorrection \
# I=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned.txt.gz \
# O=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellcorrected.txt.gz \
# CELL_BC_LIST=$my_dir/Alignments/SEEneg-cDNA/STAR_SEEneg-cDNA_readcounts.txt.gz \
# READCOUNT_CUTOFF=100 \
# SUMMARY=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellcorrected_summary.txt \
# HEADER=TRUE
#
# # UMI merging
# java -jar $my_dir/PLA_alignment.jar UMIMerging \
# I=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellcorrected.txt.gz \
# O=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellUMIcorrected.txt.gz \
# SUMMARY=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellUMIcorrected_summary.txt
#
# # Histogram of readcount (to draw kneeplot)
# java -jar $my_dir/PLA_alignment.jar ReadcountHistogram \
# I=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellcorrected.txt.gz \
# O=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellcorrected_readcounts.txt.gz

#THIS PART IN PYTHON****************************
#Get the cell barcodes from the knee plot
#Export the barcodes as text file: SEEneg-PLA_knee_barcodes.txt
#END PART IN PYTHON****************************

# Export to digital count matrix
java -jar $my_dir/PLA_alignment.jar DigitalCount \
I=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_aligned_cellUMIcorrected.txt.gz \
O=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_dge.txt.gz \
CELL_BC_LIST=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_knee_barcodes.txt \
SUMMARY=$my_dir/Alignments/SEEneg-PLA/SEEneg-PLA_dge_summary.txt
