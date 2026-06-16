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
#SBATCH -o ./out/index.%j.out             # output file
#SBATCH -e ./out/index.%j.err             # error file
#SBATCH --mail-type=END		        # Notifications (once the code words, change to END to avoid too many emails)
#SBATCH --mail-user=g.vyoshikawa@uq.edu.au

module load bowtie/1.3.1-gcc-12.3.0

cd /home/uqgventu/UMR_sorghum/genome/v5.1/assembly
gunzip Sbicolor_730_v5.0.fa.gz

bowtie-build --threads 20 Sbicolor_730_v5.0.fa
if [ $? -ne 0 ]; then
        echo "Failed to build index for Sbicolor_730_v5.0.fa"
        exit 1
    else
        echo "Index built successfully for Sbicolor_730_v5.0.fa"
    fi

gzip -v Sbicolor_730_v5.0.fa.gz