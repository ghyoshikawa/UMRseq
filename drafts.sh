#tmux session

#Log into bunya
#create new temrinal:
tmux new -s sessionName

cd /scratch/user/uqgventu/
cd /QRISdata/Q9486/data/Sb-UMRseq-diversity-panel/AGRF_CAGRF221112563_HMCLFDSX5

cp Sorghum_343_FF_SC831-14E_HMCLFDSX5_GCTTCATATT-AGGCTGAACG_L002_R*.fastq.gz /scratch/user/uqgventu/reads/
cp Sorghum_10_FF_RTx7000_HMCLFDSX5_CGGTTACGGC-CTATAGTCTT_L001_R*.fastq.gz /scratch/user/uqgventu/reads/
cp Sorghum_33_FF_SC35-14E_HMCLFDSX5_GGAATTGTAA-AGGATGTGCT_L001_R*.fastq.gz /scratch/user/uqgventu/reads/
cp Sorghum_18_R931945-2-2_HMCLFDSX5_TATTGCGCTC-CCTAACACAG_L001_R*.fastq.gz /scratch/user/uqgventu/reads/
cp Sorghum_7_FF_BTx623_HMCLFDSX5_TCCATTGCCG-TCGTGCATTC_L001_R*.fastq.gz /scratch/user/uqgventu/reads/



salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=2GB --job-name=TEST --time=02:00:00 --partition=general --qos=normal --account=a_agfs_ps srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

module load miniforge/26.1.0-0
source $ROOTMINIFORGE/etc/profile.d/conda.sh
#or 
module load miniconda3/23.9.0-0
source $EBROOTMINICONDA3/etc/profile.d/conda.sh

#check which envs are available:
conda info --envs

conda activate multiqc


cd software
pip install cutadapt  #also in multiqc env

#### lines for UMR-seq pipeline
salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=10G --job-name=TEST --time=02:00:00 --partition=general --account=a_agfs_ps srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

cd reads
find *R*.fastq.gz | sed 's/_R[12]\.fastq\.gz//' | uniq >../samples.txt

cd /scratch/user/uqgventu/
##change conda env  name with cutadapt installed
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/01-trim_galore_gz_sbatch_01.sh \
samples.txt \
20:00:00 \
16 \
cutadapt_v3.7 \
a_agfs_ps

module load miniconda3/23.9.0-0
source $EBROOTMINICONDA3/etc/profile.d/conda.sh
conda info --envs
cd /scratch/user/uqgventu/analysis/fastqc

conda activate multiqc_v3.8 || exit 1
#one line to run it
multiqc .
#after its run deactivate your environment
conda deactivate
mv multiqc_report.html multiqc_report1.html
mkdir /QRISdata/Q9486/data/multiqc_results/
mkdir /QRISdata/Q9486/data/fastqc_results/

rsync -rhiPvt multiqc_* /QRISdata/Q9486/data/multiqc_results/

cd logs/20260529-160836_01-trim_galore_gz_01
less 01-trim_galore_gz_e1
echo -e "Sample\tTotal_sequences_analysed\tPERCENT_READS_WITH_ADAPTERS_R1\tPERCENT_READS_WITH_ADAPTERS_R2\tPERCENT_BP_TRIMMED_R1\tPERCENT_BP_TRIMMED_R2" > total_reads_summary.tsv
for i in $(ls 01-trim_galore_gz_01_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
TOTAL_READS=$(grep 'Total number of sequences analysed:' $i | tr -s ' ' | cut -d " " -f 6)
PERCENT_READS_WITH_ADAPTERS=$(grep 'Reads with adapters:' $i | tr -s ' ' | cut -d " " -f 5 | paste -sd '\t')
PERCENT_BP_TRIMMED=$(grep 'Quality-trimmed:' $i | tr -s ' ' | cut -d " " -f 4 | paste -sd '\t')
echo -e "$SAMPLE\t$TOTAL_READS\t$PERCENT_READS_WITH_ADAPTERS\t$PERCENT_BP_TRIMMED"
done >> total_reads_summary.tsv

#view output:
cat total_reads_summary.tsv | column -t

cd /scratch/user/uqgventu/
vim samples_new.txt 
#OR
touch samples_new.txt
nano samples_new.txt
#make sure name order matches samples.txt to get corretn renaming in the next step
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


### QC of alignments:
# eg cd logs/20230212-202753_02-bowtie2
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


#Then calling peaks
conda create --name epic2_v3.8 python=3.8
conda activate epic2_v3.8
conda install -c bioconda epic2=0.0.41
pip install deeptools=3.5.0
deeptools --version


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


#change /analysis/trimmed_align_bowtie.
#installed pip install epic2 0.0.54
#added module load and source to eppic script


#deeptools installed with pip in poython 3.8

#add instruction to install mosdepth

#picard also doesn't work, soolve this with

module load java/17.0.15
#access using 
java -jar ~/software/picard/build/libs/picard.jar

#followed picard github instructions
#change deeptools script with 
/home/uqgventu/.conda/envs/deeptools-hacked_v3.8/lib/python3.8/site-packages/deeptools/

bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/03c-deeptools-hacked_bigWig_sbatch.sh \
samples_renamed_merged.txt \
3:00:00 \
20 \
trimmed_align_bowtie2 \
deeptools-hacked_v3.8 \
normal_res \
/home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0/Sbicolor_730_v5.0-fragments-filtered-50%Ns.bed \
a_agfs_ps

cd /home/uqgventu/gitrepos/umrseq/In-silico-digest

bash \
/home/uqgventu/gitrepos/umrseq/In-silico-digest/in-silico-digest_sbatch.sh \
-g /home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0.fa \
-o /home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0_digest \
-f 100 \
-M /home/uqgventu/UMR_sorghum/RE_motifs.fa \
-G /home/uqgventu/gitrepos/ \
-A a_agfs_ps \
-m 50 \
-c 1 \
-t 4:00:00


## to make blacklists:
bash \
/home/uqgventu/gitrepos/umrseq/In-silico-digest/in-silico-digest_sbatch.sh \
-g /home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0.fa \
-o /home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0 \
-f 100 \
-M /home/uqgventu/UMR_sorghum/RE_motifs.fa \
-G /home/uqgventu/gitrepos/umrseq/In-silico-digest \
-A a_agfs_ps \
-m 50 \
-c 1 \
-t 4:00:00 \
-e digestEnv

bash \
/home/uqgventu/gitrepos/umrseq/In-silico-digest/process-fragments_sbatch.sh \
-g /home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0.fa \
-o /home/uqgventu/UMR_sorghum/genome/v5.1/assembly/Sbicolor_730_v5.0 \
-f 100 \
-M /home/uqgventu/UMR_sorghum/RE_motifs.fa \
-G /home/uqgventu/gitrepos/ \
-A a_agfs_ps \
-m 30 \
-c 1 \
-t 10:00:00

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

#To find an error file within a logs directory for a certain sample, use:
grep -rn "your string here" /path/to/directory
grep -rn "Sorghum_233_FF_SC502-14E" .

#For FRIP
cd /scratch/user/uqgventu/analysis/trimmed_align_bowtie2_epic2_FRIP

for i in $(cat ../../samples_renamed_merged.txt); do
SAMPLE=$i
FRIP=$(awk '/percent/{getline; print}' ${i}_FRIP_counts.tab | cut -f 3)
echo -e "$SAMPLE\t$FRIP"
done > FRIP_scraped.tsv
# add headder
echo -e 'sample\tFRIP_score' \
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

sed '1d' FRIP_scraped.tsv |
awk '{ total += $2 } END { print total/NR }' -

## but Akshya says at this read depth we should be getting double the peaks.
#Akshya runs epic2 with bin size 50. Lets try that next.