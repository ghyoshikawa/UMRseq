module load miniconda3/23.9.0-0
source $EBROOTMINICONDA3/etc/profile.d/conda.sh

bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/01-trim_galore_gz_sbatch_01.sh \
samples.txt \
20:00:00 \
16 \
cutadapt_v3.7 \
a_agfs_ps

less 01-trim_galore_gz_e1
echo -e "Sample\tTotal_sequences_analysed\tPERCENT_READS_WITH_ADAPTERS_R1\tPERCENT_READS_WITH_ADAPTERS_R2\tPERCENT_BP_TRIMMED_R1\tPERCENT_BP_TRIMMED_R2" > total_reads_summary.tsv
for i in $(ls 01-trim_galore_gz_01_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
TOTAL_READS=$(grep 'Total number of sequences analysed:' $i | tr -s ' ' | cut -d " " -f 6)
PERCENT_READS_WITH_ADAPTERS=$(grep 'Reads with adapters:' $i | tr -s ' ' | cut -d " " -f 5 | paste -sd '\t')
PERCENT_BP_TRIMMED=$(grep 'Quality-trimmed:' $i | tr -s ' ' | cut -d " " -f 4 | paste -sd '\t')
echo -e "$SAMPLE\t$TOTAL_READS\t$PERCENT_READS_WITH_ADAPTERS\t$PERCENT_BP_TRIMMED"
done >> total_reads_summary.tsv
cat total_reads_summary.tsv | column -t


sort samples.txt > samples_sorted.txt
sort samples_new.txt > samples_new_sorted.txt

paste -d '\t' samples_sorted.txt <(cut -f1 samples_new_sorted.txt) > samples_key.tsv
cat samples_key.tsv | column -t
mv analysis/trimmed analysis/trimmed_lanes
mkdir analysis/trimmed
sbatch --array=1-$(wc -l < samples_key.tsv) cat.sh samples_key.tsv

#alignment:
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/02-bowtie2_sbatch_01.sh \
samples_renamed_merged.txt \
/scratch/user/uqgventu/analysis/trimmed/ \
6 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0 \
10 \
18:00:00 \
40 \
a_agfs_ps

cd /scratch/user/uqgventu/logs/20260617-175654_02-bowtie2

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


bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/05-epic2_sbatch.sh \
samples_renamed_merged.txt \
03:00:00 \
20 \
/scratch/user/uqgventu/analysis/trimmed_align_bowtie2 \
epic2_v3.8 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_730_v5.0.chrom.sizes \
100 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_730_v5.0-fragments-filtered-50%Ns.bed \
a_agfs_ps

#To get summary data from epic2 output, run the following command:
bin_Size=100

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


# should look like this (except my fragments are pretty small...)
sample                     MEDIAN_LENGTH  MEAN_LENGTH  NUM_REGIONS
Sorghum_392_FF_SC1063-14E  1199           1754.82      83939
Sorghum_441_FF_SC1260-5
Sorghum_153_FF_SC325-14E   899            1246.21      97342
Sorghum_249_FF_SC556-3     999            1447.37      92702
Sorghum_463_FF_SC1325-14E  1299           1719.79      70103

# and here is some code to get the averages:
sed '1d' epic2_region_metrics.tsv |
awk '{ total += $2 } END { print total/NR }' -

# median frag size
1113.14

sed '1d' epic2_region_metrics.tsv |
awk '{ total += $3 } END { print total/NR }' -

# mean frag size
1550.09
#mean nr of regions
sed '1d' epic2_${bin_Size}_region_metrics.tsv |
awk '{ total += $4 } END { print total/NR }' -
#80959.7

#To find an error file within a logs directory for a certain sample, use:
grep -rn "your string here" /path/to/directory
grep -rn "Sorghum_233_FF_SC502-14E" .

#OOM
: '
Sorghum_233_FF_SC502-14E.sam
123936632 reads; of these:
  123936632 (100.00%) were paired; of these:
    6058635 (4.89%) aligned concordantly 0 times
    61565765 (49.68%) aligned concordantly exactly 1 time
    56312232 (45.44%) aligned concordantly >1 times
    ----
    6058635 pairs aligned concordantly 0 times; of these:
      466528 (7.70%) aligned discordantly 1 time
    ----
    5592107 pairs aligned 0 times concordantly or discordantly; of these:
      11184214 mates make up the pairs; of these:
        8861018 (79.23%) aligned 0 times
        832443 (7.44%) aligned exactly 1 time
        1490753 (13.33%) aligned >1 times
96.43% overall alignment rate
+ echo 'Sorghum_233_FF_SC502-14E total alignments before MAPQ filter (might include reads that map to multiple locations)'
+ samtools view -c /scratch/user/uqgventu/analysis/trimmed/_align_bowtie2/Sorghum_233_FF_SC502-14E.sam
+ samtools view -q 10 -b -@ 6 /scratch/user/uqgventu/analysis/trimmed/_align_bowtie2/Sorghum_233_FF_SC502-14E.sam
'
#For FRIP
cd /scratch/user/uqgventu/analysis/trimmed_align_bowtie2_epic2_FRIP

for i in $(cat ../../samples_renamed_merged.txt); do
SAMPLE=$i
FRIP=$(awk '/percent/{getline; print}' ${i}_FRIP_counts.tab | cut -f 3)
echo -e "$SAMPLE\t$FRIP"
done > FRIP_scraped.tsv
# add headder
echo -e 'sample\tFRIP score' \
| cat - FRIP_scraped.tsv > temp && mv temp FRIP_scraped.tsv

cat FRIP_scraped.tsv | column -t

#Looks ok
sample                     FRIP    score
Sorghum_392_FF_SC1063-14E  89.42
Sorghum_441_FF_SC1260-5
Sorghum_153_FF_SC325-14E   86.74
Sorghum_249_FF_SC556-3     90.86
Sorghum_463_FF_SC1325-14E  88.85
Sorghum_314_FF_SC715-14E   89.92
Sorghum_312_FF_SC708-14E   91.10

## but Akshya says at this read depth we should be getting double the peaks.
#Akshya runs epic2 with bin size 50. Lets try that next.
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/05-epic2_sbatch.sh \
samples_renamed_merged.txt \
03:00:00 \
20 \
/scratch/user/uqgventu/analysis/trimmed_align_bowtie2 \
epic2_v3.8 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_730_v5.0.chrom.sizes \
50 \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_730_v5.0-fragments-filtered-50%Ns.bed \
a_agfs_ps

cd /scratch/user/uqgventu/analysis/trimmed_align_bowtie2_epic2
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
done > epic2_${bin_Size}_region_metrics.tsv

# add header
echo -e 'sample\tMEDIAN_LENGTH\tMEAN_LENGTH\tNUM_REGIONS' \
| cat - epic2_${bin_Size}_region_metrics.tsv > temp && mv temp epic2_${bin_Size}_region_metrics.tsv

cat epic2_${bin_Size}_region_metrics.tsv | column -t

# and here is some code to get the median frag size:
sed '1d' epic2_${bin_Size}_region_metrics.tsv |
awk '{ total += $2 } END { print total/NR }' -
#659.7
# mean frag size
sed '1d' epic2_${bin_Size}_region_metrics.tsv |
awk '{ total += $3 } END { print total/NR }' -
#884.9
#number of regions
sed '1d' epic2_${bin_Size}_region_metrics.tsv |
awk '{ total += $4 } END { print total/NR }' -
#117620

#Reducing the bin size reduced the median fragment size from 1113.14 to 659.7 
#and the mean fragment size from 1550.09 to 884.9, nr of regions from 80959.7 to 117620. 
#So this is a good thing, better resolution

#FRIP score with bin size 50 (bin size 100 got overwritten)

for i in $(cat ../../samples_renamed_merged.txt); do
SAMPLE=$i
FRIP=$(awk '/percent/{getline; print}' ${i}_FRIP_counts.tab | cut -f 3)
echo -e "$SAMPLE\t$FRIP"
done > FRIP_scraped.tsv
# add headder
echo -e 'sample\tFRIP score' \
| cat - FRIP_scraped.tsv > temp && mv temp FRIP_scraped.tsv

cat FRIP_scraped.tsv | column -t
sed '1d' FRIP_scraped.tsv |
awk '{ total += $2 } END { print total/NR }' -
#average FRIP score 84.37