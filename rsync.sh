#!/bin/bash --login
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1            #change for multi-threaded jobs
#SBATCH --array=1-1                 #change for array jobs, e.g. 1-10 for 10 tasks  
#SBATCH --time=1-00:60:00             # time allocation (D-HH:MM:SS)
#SBATCH --mem=20G                   # memory pool for all cores
#SBATCH --qos=normal
#SBATCH --account=a_agfs_ps
#SBATCH --job-name=rsync
#SBATCH -o ./out/rsync.%j.out             # output file
#SBATCH -e ./out/rsync.%j.err             # error file
#SBATCH --mail-type=END		        # Notifications (once the code words, change to END to avoid too many emails)
#SBATCH --mail-user=g.vyoshikawa@uq.edu.au

rsync AGRF_CAGRF221112563_HMCLFDSX5.tar /QRISdata/Q9486/data/
