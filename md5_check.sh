#!/bin/bash --login
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1            #change for multi-threaded jobs
#SBATCH --array=1-20                 #change for array jobs, e.g. 1-10 for 10 tasks  
#SBATCH --time=1-00:00:00             # time allocation (D-HH:MM:SS)
#SBATCH --mem=100G                   # memory pool for all cores
#SBATCH --qos=normal
#SBATCH --account=a_agfs_ps
#SBATCH --job-name=CHANGE-ME
#SBATCH -o ./out/md5.%j.out             # output file
#SBATCH -e ./out/md5.%j.err             # error file
#SBATCH --mail-type=END		        # Notifications (once the code words, change to END to avoid too many emails)
#SBATCH --mail-user=g.vyoshikawa@uq.edu.au

module purge
module load parallel/20230722-gcccore-12.3.0

DIR=/QRISdata/Q9486/data/Sb-UMRseq-diversity-panel/AGRF_CAGRF221112563_HMCLFDSX5/
cd $DIR
# Create new output file with timestamp

TIMESTAMP=$(date +"%Y%m%d")
OUTPUT_FILE="md5hashes_merged_${TIMESTAMP}.txt"
echo "Creating MD5 hashes in: $OUTPUT_FILE"
touch "$OUTPUT_FILE"

export OUTPUT_FILE
export DIR

md5_merged() {
    SAMPLE="$1"
md5sum -c checksums.md5

if [ $? -ne 0 ]; then
        echo "Failed to calculate MD5 for $DIR/$1"
        exit 1
    else
        echo "MD5 calculated successfully for $DIR/$1"
    fi
}
export -f md5_merged  # Export so GNU parallel can use it

ls *.fastq.gz | sort | uniq | parallel -j 20 md5_merged {}
if [ $? -ne 0 ]; then
    echo "Failed to list files in $DIR"
    exit 1
fi