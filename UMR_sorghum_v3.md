# running the UMRseq pipeline for sorghum population data
# using previous trimmed reads, to align to v3 for Dr Claire Huang

# for v3.0.1 reference genome:
github repo has been recloned to account for previous changes. File names have previously been edited to group lanes and drop adapter sequences from file name.

```bash
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

cd /scratch/user/uqgventu/logs/20261026-143823_02-bowtie2

# scrape logs to get summary and also to check mapping of all samples worked
echo -e "sample\tALIGNED_1_TIME\tMULTI_MAPPINGS\tUNMAPPED" > bowtie2_summary.tsv

for i in $(ls 02-bowtie2_e*); do
SAMPLE=$(grep 'echo sample being mapped is' $i | cut -d " " -f 7)
ALIGNED_1_TIME=$(grep ') aligned concordantly exactly 1 time' $i | cut -d " " -f 6)
MULTI_MAPPINGS=$(grep ' aligned concordantly >1 times' $i | cut -d " " -f 6)
UNMAPPED=$(grep ') aligned concordantly 0 times' $i | cut -d " " -f 6)
echo -e "$SAMPLE\t$ALIGNED_1_TIME\t$MULTI_MAPPINGS\t$UNMAPPED"
done >> bowtie2_summary.tsv

cat bowtie2_summary.tsv | column -t 

echo -e "sample\tTOTAL_ALIGNMENTS\tMAPQ10\tPERCENT" >MAPQ_filter_summary.tsv

for i in $(ls 02-bowtie2_o*); do
SAMPLE=$(grep 'sample being mapped is' $i | cut -d " " -f 5)
TOTAL_ALIGNMENTS=$(grep -A1 'total alignments before MAPQ filter' $i | grep -v "total alignments before MAPQ filter")
MAPQ10=$(grep -A1 'total alignments after MAPQ filter' $i | grep -v "total alignments after MAPQ filter")
PERCENT=$(bc -l <<< "$MAPQ10/$TOTAL_ALIGNMENTS*100" | awk '{printf("%d\n",$1 + 0.5)}')
echo -e "$SAMPLE\t$TOTAL_ALIGNMENTS\t$MAPQ10\t$PERCENT"
done >> MAPQ_filter_summary.tsv

cat MAPQ_filter_summary.tsv | column -t 
```

# Have to create a blacklist for v3 oops

```bash 
bash \
/home/$USER/gitrepos/umrseq/In-silico-digest/in-silico-digest_sbatch.sh \
-g /home/uqgventu/UMR_sorghum/genome/v3.0/assembly/Sbicolor_454_v3.0.1.fa \
-o /home/uqgventu/UMR_sorghum/genome/v3.0/assembly/output \
-f 100 \
-M /home/$USER/UMR_sorghum/RE_motifs.fa \
-G /home/$USER/gitrepos/umrseq/In-silico-digest \
-A a_agfs_ps \
-m 50 \
-c 1 \
-t 4:00:00 \
-e digestEnv

```

# Now call umrs using epic2, change to v3 files

```bash
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/05-epic2_sbatch.sh \
samples_renamed_merged.txt \
03:00:00 \
20 \
/scratch/user/uqgventu/analysis/trimmed_align_bowtie2 \
epic2_v3.8 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_454_v3.0.1.chrom.sizes \
100 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_730_v5.0-fragments-filtered-50%Ns.bed \
a_agfs_ps

