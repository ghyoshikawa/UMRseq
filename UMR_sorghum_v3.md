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
cd /home/uqgventu/UMR_sorghum/genome/v3.0/assembly/

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

bash \
/home/uqgventu/gitrepos/umrseq/In-silico-digest/process-fragments_sbatch.sh \
-g /home/uqgventu/UMR_sorghum/genome/v3.0/assembly/Sbicolor_454_v3.0.1.fa \
-o /home/uqgventu/UMR_sorghum/genome/v3.0/assembly/output \
-f 100 \
-M /home/uqgventu/UMR_sorghum/RE_motifs.fa \
-G /home/$USER/gitrepos/umrseq/In-silico-digest \
-A a_agfs_ps \
-m 30 \
-c 1 \
-t 15:00:00
#this block failed with 10h walltime, seems like it only had a few more lines to go. increased walltime to 15h

```

# Now call umrs using epic2, change to v3 files

```bash
cd /scratch/user/uqgventu
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/05-epic2_sbatch.sh \
samples_renamed_merged.txt \
03:00:00 \
20 \
/scratch/user/uqgventu/analysis/trimmed_align_bowtie2 \
epic2_v3.8 \
/home/uqgventu/UMR_sorghum/genome/v3.0/assembly/output/Sbicolor_454_v3.0.1.chrom.sizes \
50 \
/home/uqgventu/UMR_sorghum/genome/v3.0/assembly/output/Sbicolor_454_v3.0.1-fragments-filtered-50%Ns.bed \
a_agfs_ps

bin_Size=50

# set this to match your epic2 output naming
for i in $(cat ../../samples_renamed_merged.txt); do
  SAMPLE=$i
  BEDFILE="${i}.${bin_Size}.epic2.bed"

  METRICS=$(awk '
    {
      len = $3 - $2          # region length = End - Start
      sum += len
      lengths[NR] = len
      n++
    }
    END {
      if (n == 0) { print "NA\tNA\t0"; exit }   # guard: no regions
      # median
      asort(lengths)
      if (n % 2) median = lengths[(n+1)/2]
      else       median = (lengths[n/2] + lengths[n/2+1]) / 2
      # mean, count
      printf "%s\t%s\t%s", median, sum/n, n
    }' "$BEDFILE")

  echo -e "$SAMPLE\t$METRICS"
done > epic2_region_metrics.tsv

# add header
echo -e 'sample\tMEDIAN_LENGTH\tMEAN_LENGTH\tNUM_REGIONS' \
| cat - epic2_region_metrics.tsv > temp && mv temp epic2_region_metrics.tsv

cat epic2_region_metrics.tsv | column -t 

# and here is some code to get the averages:
#median
sed '1d' epic2_region_metrics.tsv |
awk '{ total += $2 } END { print total/NR }' -
#mean
sed '1d' epic2_region_metrics.tsv |
awk '{ total += $3 } END { print total/NR }' -

```

# consensus UMRs 
remove low quality samples
```bash
#remove low-read samples from epic2 analysis, FRiP <60, >9M reads

tail -n +2 low_qual_samples.csv | cut -d',' -f1 | \
grep -vFf - samples_renamed_merged.txt > samples_renamed_merged_filtered.txt

cd /scratch/user/$USER/

#change the file suffix accordingly
bash ~/gitrepos/umrseq/consensus_UMR/20-EPIC2_reformat_sbatch.sh \
samples_renamed_merged_filtered.txt \
.50.epic2.bed \
./analysis/trimmed_align_bowtie2_epic2 \
./analysis/umr_names \
20:00 \
2 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/21-UMR_consensus_merge_sbatch.sh  \
samples_renamed_merged_filtered.txt \
.50.epic2.bed \
./analysis/umr_names \
./analysis/umr_consensus \
all_samples \
50 \
30:00 \
20 \
UMR_calling \
a_agfs_ps