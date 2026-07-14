#!/bin/bash --login
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1            #change for multi-threaded jobs
#SBATCH --time=2-00:00:00             # time allocation (D-HH:MM:SS)
#SBATCH --mem=10                   # memory pool for all cores
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

md5sum -c --ignore-missing checksums.md5

if [ $? -ne 0 ]; then
    echo "MD5 check failed for one or more files. Please investigate."
    exit 1
else
    echo "All files passed the MD5 check successfully."
fi