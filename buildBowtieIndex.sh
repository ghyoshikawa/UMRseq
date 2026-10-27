#!/bin/bash --login
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20            #change for multi-threaded jobs
#SBATCH --time=1-00:00:00             # time allocation (D-HH:MM:SS)
#SBATCH --mem=100G                   # memory pool for all cores
#SBATCH --qos=normal
#SBATCH --account=a_agfs_ps
#SBATCH --job-name=buildBowtieIndex           # sensible name for the job
#SBATCH -o ./out/index.%j.out             # output file
#SBATCH -e ./out/index.%j.err             # error file
#SBATCH --mail-type=END		        # Notifications (once the code words, change to END to avoid too many emails)
#SBATCH --mail-user=g.vyoshikawa@uq.edu.au

module load bowtie2/2.5.1-gcc-12.3.0 samtools/1.18-gcc-12.3.0

cd /home/uqgventu/UMR_sorghum/genome/v5.1/assembly
 if [ -f Sbicolor_730_v5.0.fa.gz ]; then
   echo "File Sbicolor_730_v5.0.fa.gz found, uncompressing..."
   gunzip Sbicolor_730_v5.0.fa.gz
   elif [ -f Sbicolor_730_v5.0.fa ]; then
     echo "File Sbicolor_730_v5.0.fa found and uncompressed successfully"
   else
     echo "No genome found"
fi

#use --large-index if the genome is larger than 4GB
bowtie2-build --threads 20 Sbicolor_730_v5.0.fa Sbicolor_730_v5.0
if [ $? -ne 0 ]; then
        echo "Failed to build index for Sbicolor_730_v5.0.fa"
        exit 1
    else
        echo "Index built successfully for Sbicolor_730_v5.0.fa"
    fi

bowtie2-inspect -vs Sbicolor_730_v5.0 > Sbicolor_730_v5.0_index_info.txt

#Define chromosome sizes
samtools faidx Sbicolor_730_v5.0.fa
cut -f1,2 Sbicolor_730_v5.0.fa.fai > Sbicolor_730_v5.0_chrom.sizes

gzip -v Sbicolor_730_v5.0.fa