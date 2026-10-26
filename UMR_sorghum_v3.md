# running the UMRseq pipeline for sorghum population data
# using previous trimmed reads, to align to v3 for Dr Claire Huang

# for v3.0.1 reference genome:
github repo has been recloned to account for previous changes. File names have previously been edited to group lanes and drop adapter sequences from file name.

```
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/02-bowtie2_sbatch.sh \
samples_renamed_merged.txt \
/scratch/user/uqgventu/analysis/trimmed/ \
6 \
~/refseqs/sorghum/Sbicolor_454_v3.0.1 \
10 \
18:00:00 \
40 \
a_agfs_ps

```

