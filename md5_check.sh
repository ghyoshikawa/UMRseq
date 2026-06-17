#!/bin/bash --login
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20            #change for multi-threaded jobs
#SBATCH --time=1-00:00:00             # time allocation (D-HH:MM:SS)
#SBATCH --mem=50G                   # memory pool for all cores
#SBATCH --qos=normal
#SBATCH --account=a_agfs_ps
#SBATCH --job-name=md5_check           # sensible name for the job
#SBATCH -o ./out/md5.%j.out             # output file
#SBATCH -e ./out/md5.%j.err             # error file
#SBATCH --mail-type=END		        # Notifications (once the code words, change to END to avoid too many emails)
#SBATCH --mail-user=g.vyoshikawa@uq.edu.au

module purge
module load parallel/20230722-gcccore-12.3.0

DIR=/scratch/user/uqgventu/AGRF_CAGRF221112563_HMCLFDSX5/
cd $DIR
# Create new output file with timestamp

TIMESTAMP=$(date +"%Y%m%d")
export DIR

md5_merged() {
    SAMPLE="$1"
md5sum -c --ignore-missing checksums.md5

}
export -f md5_merged  # Export so GNU parallel can use it

ls *.fastq.gz | sort | uniq | parallel -j 20 md5_merged {}
if [ $? -ne 0 ]; then
    echo "Failed to list files in $DIR"
    exit 1
fi