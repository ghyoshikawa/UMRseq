# running the UMRseq pipeline for sorghum population data
# using previous trimmed reads, to align to v3 for Dr Claire Huang

For v3.0.1 reference genome:
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
Skip this bit, I gave up as I found a file 

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
-t 1-0:00:00
#this block failed with 10/15h walltime. 


```

Now call umrs using epic2, change to v3 files, bin size 50 bp

```bash
cd /scratch/user/uqgventu
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/05-epic2_sbatch.sh \
samples_renamed_merged.txt \
03:00:00 \
20 \
/scratch/user/uqgventu/analysis/trimmed_align_bowtie2 \
epic2_v3.8 \
~/UMR_sorghum/genome/v3.0/assembly/output/Sbicolor_454_v3.0.1.chrom.sizes \
50 \
~/UMR_sorghum/genome/v3.0/assembly/output/Sbicolor_454_v3.0.1-fragments-filtered-50%Ns.bed \
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
#650 bp
#mean
sed '1d' epic2_region_metrics.tsv |
awk '{ total += $3 } END { print total/NR }' -
#870 bp

#scrape FRIP logs:
cd /scratch/user/$USER/analysis/trimmed_align_bowtie2_epic2_FRIP
for i in $(cat ../../samples_renamed_merged.txt); do
  SAMPLE=$i
  FRIP=$(awk '/percent/{getline; print}' ${i}_FRIP_counts.tab | cut -f 3)
  echo -e "$SAMPLE\t$FRIP"
  done > FRIP_scraped.tsv
# add headder
echo -e 'sample\tFRIP score' \
| cat - FRIP_scraped.tsv > temp && mv temp FRIP_scraped.tsv

cat FRIP_scraped.tsv | column -t

#get FRIP score average:
sed '1d' FRIP_scraped.tsv |
awk '{ total += $2 } END { print total/NR }' -
```

# consensus UMRs for v3
remove low quality samples
```bash
#remove low-read samples from epic2 analysis, FRiP <60, >9M reads
cd ../..

tail -n +2 low_qual_v3.csv | cut -d',' -f1 | \
grep -vFf - samples_renamed_merged.txt > samples_renamed_merged_filtered.txt

#25 samples lost due to <9M reads or FRiP<60

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
100 \
30:00 \
20 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/21b-UMR_consensus_merge_background_sbatch.sh \
./analysis/umr_consensus/all_samples.umr.merged.consensus.noALT.bed \
~/UMR_sorghum/genome/v3.0/assembly/Sbicolor_454_v3.0.1.fa.fai \
./analysis/umr_consensus \
0.2 \
12345 \
2:00:00 \
8 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/22-UMR_featurecounts_sbatch.sh \
samples_renamed_merged_filtered.txt \
.bam \
./analysis/umr_consensus/all_samples.umr.merged.consensus.noALT.saf \
./analysis/trimmed_align_bowtie2 \
./analysis/umr_counts \
20:00 \
6 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/23-UMR_featurecounts_summary_sbatch.sh \
samples_renamed_merged_filtered.txt \
./analysis/umr_counts \
./analysis/umr_counts \
all_samples \
0:15:00 \
2 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/24-UMR_combine_calls_sbatch.sh \
samples_renamed_merged_filtered.txt \
./analysis/umr_counts \
./analysis/umr_calling \
all_samples \
1000 \
15:00 \
5 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/25-UMR_calls_sbatch.sh \
./analysis/umr_calling/batches/batch_list.txt \
./analysis/umr_calling/calls \
30:00 \
5 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/26-UMR_calls_table_sbatch.sh \
./analysis/umr_calling/batches/batch_list.txt \
samples_renamed_merged_filtered.txt \
./analysis/umr_calling/calls \
./analysis/umr_calling \
all_samples \
0:10:00 \
2 \
a_agfs_ps

#Done!
```
# consensus UMRs for v5 with epic2 bin size 100bp

remove low quality samples
```bash
#remove low-read samples from epic2 analysis, FRiP <60, >9M reads

tail -n +2 low_qual_samples.csv | cut -d',' -f1 | \
grep -vFf - samples_renamed_merged.txt > samples_renamed_merged_filtered.txt

cd /scratch/user/$USER/v5

#change the file suffix accordingly
bash ~/gitrepos/umrseq/consensus_UMR/20-EPIC2_reformat_sbatch.sh \
samples_renamed_merged_filtered.txt \
.100.epic2.bed \
./analysis/trimmed_align_bowtie2_epic2_bin100 \
./analysis/umr_names \
20:00 \
2 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/21-UMR_consensus_merge_sbatch.sh  \
samples_renamed_merged_filtered.txt \
.100.epic2.bed \
./analysis/umr_names \
./analysis/umr_consensus \
all_samples \
100 \
30:00 \
20 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/21b-UMR_consensus_merge_background_sbatch.sh \
./analysis/umr_consensus/all_samples.umr.merged.consensus.noALT.bed \
~/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0.fa.fai \
./analysis/umr_consensus \
0.2 \
12345 \
2:00:00 \
8 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/22-UMR_featurecounts_sbatch.sh \
samples_renamed_merged_filtered.txt \
.bam \
./analysis/umr_consensus/all_samples.umr.merged.consensus.noALT.saf \
./analysis/trimmed_align_bowtie2 \
./analysis/umr_counts \
20:00 \
6 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/23-UMR_featurecounts_summary_sbatch.sh \
samples_renamed_merged_filtered.txt \
./analysis/umr_counts \
./analysis/umr_counts \
all_samples \
0:15:00 \
2 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/24-UMR_combine_calls_sbatch.sh \
samples_renamed_merged_filtered.txt \
./analysis/umr_counts \
./analysis/umr_calling \
all_samples \
1000 \
15:00 \
5 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/25-UMR_calls_sbatch.sh \
./analysis/umr_calling/batches/batch_list.txt \
./analysis/umr_calling/calls \
30:00 \
5 \
UMR_calling \
a_agfs_ps

bash ~/gitrepos/umrseq/consensus_UMR/26-UMR_calls_table_sbatch.sh \
./analysis/umr_calling/batches/batch_list.txt \
samples_renamed_merged_filtered.txt \
./analysis/umr_calling/calls \
./analysis/umr_calling \
all_samples \
0:10:00 \
2 \
a_agfs_ps

```
# SNP investigation from UMRs

Now that the alignments and consensus UMRs have been called, we proceed to call SNPs from UMRs. 

Also have to look into generating consensus UMRs based on read count levels rather than UMR presence or absence

Both Pete and Claire have shared a script, so can compare both of them and pick the best?


```