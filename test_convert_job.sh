#!/bin/bash
#SBATCH --job-name AmpliconTest3_ConvertJob
#SBATCH -p qiita
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time 36:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user c2cowart@ucsd.edu
#SBATCH --mem-per-cpu 10gb
set -x
date
hostname
cd /sequencing/seqmount/KL_MiSeq_Runs/200226_M05314_0237_000000000-CMNW8
module load bclconvert_3.7.5
bcl-convert --sample-sheet "/pscratch/seq_test/AmpliconDevelopment/test_sample_sheet.csv" --output-directory "/pscratch/seq_test/AmpliconDevelopment/ConvertJob" --bcl-input-directory . --bcl-num-decompression-threads 8 --bcl-num-conversion-threads 8 --bcl-num-compression-threads 16 --bcl-num-parallel-tiles 8 --bcl-sampleproject-subdirectories true --force

