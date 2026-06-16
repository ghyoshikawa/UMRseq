#!/bin/bash --login
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1            #change for multi-threaded jobs
#SBATCH --array=1-20                 #change for array jobs, e.g. 1-10 for 10 tasks  
#SBATCH --time=00:60:00             # time allocation (D-HH:MM:SS)
#SBATCH --mem=10G                   # memory pool for all cores
#SBATCH --qos=normal
#SBATCH --account=a_agfs_ps
#SBATCH --job-name=CHANGE-ME
#SBATCH -o ./out/cat.%j.out             # output file
#SBATCH -e ./out/cat.%j.err             # error file
#SBATCH --mail-type=END		        # Notifications (once the code words, change to END to avoid too many emails)
#SBATCH --mail-user=g.vyoshikawa@uq.edu.au

# Step 1: Backup and prepare
mv analysis/trimmed analysis/trimmed_lanes
mkdir analysis/trimmed

# Step 2: Merge FASTQ files
while read -r i; do
    orgfilename=$(echo "$i" | awk '{print $1}')
    newfilename=$(echo "$i" | awk '{print $2}')
    
    # Remove the lane suffix (_L001, _L002, etc.) from orgfilename
    # This leaves just the base sample name without lane info
    basefilename=$(echo "$orgfilename" | sed 's/_L[0-9][0-9][0-9]$//')
    
    echo "Merging: $basefilename (all lanes) → $newfilename"
    
    # Now wildcard matches ALL lanes: _L001, _L002, _L003, etc.
    cat "analysis/trimmed_lanes/${basefilename}"_L*_R1_val_1.fq.gz \
        > "analysis/trimmed/${newfilename}_R1_val_1.fq.gz"
    
    cat "analysis/trimmed_lanes/${basefilename}"_L*_R2_val_2.fq.gz \
        > "analysis/trimmed/${newfilename}_R2_val_2.fq.gz"
    
done < samples_key.tsv

# Step 3: Create new sample list
cd analysis/trimmed || exit 1
find *R1_val_1.fq.gz | \
    sed 's|_R1_val_1\.fq\.gz||' \
    > ../../samples_rename_merge.txt
cd ../..

# Step 4: Back up old sample list
mv samples.txt samples_old.txt

echo "Done!"