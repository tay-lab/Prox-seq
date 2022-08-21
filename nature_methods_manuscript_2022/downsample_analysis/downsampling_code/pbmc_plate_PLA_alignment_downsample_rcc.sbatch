#!/bin/bash
#SBATCH --job-name=SmartPLA_alignment_downsample
#SBATCH --output=/project2/tays/Van/20201202_smart_PLA/sbatch_files/SmartPLA_alignment_downsample_sbatch.out
#SBATCH --error=/project2/tays/Van/20201202_smart_PLA/sbatch_files/SmartPLA_alignment_downsample_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G

#Directory variable of working Directory
my_dir=/project2/tays/Van/20201202_smart_PLA

# Directory for storing summary files

mkdir -p "$my_dir"/FastQ_downsample

# Download seqkit for read downsampling
# wget https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz
seqkit=/project2/tays/Van/seqkit/seqkit
rm -f "$my_dir"/downsample_seqkit.out
rm -f "$my_dir"/downsample_seqkit.err

# Load java**
module load java

for i in 0.005 0.01 0.05 0.1 0.25 0.5 0.75; do
  # Reformat i: replace . with _
  j="${i/./_}"
  mkdir -p "$my_dir"/downsample_"$j"

  # Record number of reads
  echo "name,reads" > "$my_dir"/downsample_"$j".csv
  echo $(date) ===== "$i" >> "$my_dir"/downsample_seqkit.out

  # Downsample
  for k in "$my_dir"/FastQ/*.fastq.gz; do
    temp1="${k/FastQ/FastQ_downsample}"
    $seqkit sample $k \
      -p "$i" -s 0 -o "$temp1" \
      >> /dev/null 2>> "$my_dir"/downsample_seqkit.err

    # Count reads after downsampling
    temp2=$(zcat "$temp1" | wc -l)
    echo "$temp1",$(($temp2 / 4)) >> "$my_dir"/downsample_"$j".csv
  done
  echo "" >> "$my_dir"/downsample_seqkit.err

  # Align raw reads
  java -jar "$my_dir"/PLA_alignment.jar ReadAlignmentSmartSeq \
  R1_LIST="$my_dir"/PBMC_R1List_downsample.csv \
  O="$my_dir"/downsample_"$j"/ReadAlignmentSmartSeq_out.txt.gz \
  AB_BC_LIST="$my_dir"/barcode_cocktail_2.6_alignment.csv \
  SUMMARY="$my_dir"/downsample_"$j"/ReadAlignmentSmartSeq_summary.txt \
  HEADER=TRUE \
  BASE_QUALITY=10 \
  NUM_BELOW_BASE_QUALITY=1 \
  BASE_QUALITY_START=1

  # UMI merging
  java -jar $my_dir/PLA_alignment.jar UMIMerging \
  I=$my_dir/downsample_"$j"/ReadAlignmentSmartSeq_out.txt.gz \
  O=$my_dir/downsample_"$j"/UMIMerging_out.txt.gz \
  SUMMARY=$my_dir/downsample_"$j"/UMIMerging_summary.txt

  # Digital count matrix
  java -jar "$my_dir"/PLA_alignment.jar DigitalCount \
  I="$my_dir"/downsample_"$j"/UMIMerging_out.txt.gz \
  O="$my_dir"/downsample_"$j"/PLA_count_matrix.txt.gz \
  CELL_BC_LIST=NONE \
  DUPLICATE_EXPORT="$my_dir"/downsample_"$j"/DigitalCountduplicate_export.txt.gz \
  REMOVE_DUPLICATE=TRUE \
  SUMMARY="$my_dir"/downsample_"$j"/DigitalCount_summary.txt

  # Delete old downsample files
  rm "$my_dir"/FastQ_downsample/*.fastq.gz
done
