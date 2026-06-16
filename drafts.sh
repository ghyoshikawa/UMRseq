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



salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=10G --job-name=TEST --time=02:00:00 --partition=general --qos=normal --account=a_agfs_ps srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l

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

conda activate multiqc_v3.8 
#one line to run it
multiqc .
#after its run deactivate your environment
conda deactivate
mv multiqc_report.html multiqc_report1.html
rsync -rhiPvt multiqc_* /QRISdata/Q9486/data/

cd logs/20260529-160836_01-trim_galore_gz_01
less 01-trim_galore_gz_e1
echo -e "Sample\tTotal_sequences_analysed\tPERCENT_READS_WITH_ADAPTERS_R1\tPERCENT_READS_WITH_ADAPTERS_R2\tPERCENT_BP_TRIMMED_R1\tPERCENT_BP_TRIMMED_R2" > ../total_reads_summary1.tsv
for i in $(ls 01-trim_galore_gz_01_e*); do
SAMPLE=$(grep '+ ID=' $i | cut -d "=" -f 2)
TOTAL_READS=$(grep 'Total number of sequences analysed:' $i | tr -s ' ' | cut -d " " -f 6)
PERCENT_READS_WITH_ADAPTERS=$(grep 'Reads with adapters:' $i | tr -s ' ' | cut -d " " -f 5 | paste -sd '\t')
PERCENT_BP_TRIMMED=$(grep 'Quality-trimmed:' $i | tr -s ' ' | cut -d " " -f 4 | paste -sd '\t')
echo -e "$SAMPLE\t$TOTAL_READS\t$PERCENT_READS_WITH_ADAPTERS\t$PERCENT_BP_TRIMMED"
done >> ../total_reads_summary1.tsv

#view output:
cat ../total_reads_summary1.tsv | column -t

cd /scratch/user/uqgventu/
vim samples_new.txt 
#OR
touch samples_new.txt
nano samples_new.txt
#make sure name order matches samples.txt to get corretn renaming in the next step

paste -d '\t' samples.txt <(cut -f1 samples_new.txt) > samples_key.tsv
cat samples_key.tsv | column -t
mv analysis/trimmed analysis/trimmed_lanes
mkdir analysis/trimmed
sbatch cat.sh

#alignment:
bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/02-bowtie2_sbatch_01.sh \
samples_rename_merge.txt \
trimmed \
6 \
~/refseqs/sorghum/Sbicolor_454_v3.0.1 \
10 \
18:00:00 \
40 \
a_agfs_ps


### QC of alignments:

#Then calling peaks
conda create --name epic2_v3.8 python=3.8
conda activate epic2_v3.8
conda install -c bioconda epic2=0.0.41
pip install deeptools=3.5.0
deeptools --version


bash \
/home/uqgventu/gitrepos/umrseq/UMRseq/05-epic2_sbatch.sh \
samples_rename_merge.txt \
0:25:00 \
20 \
/scratch/user/uqgventu/analysis/trimmed_align_bowtie2 \
epic2_v3.8 \
/home/uqgventu/refseqs/sorghum/Sbicolor_454_v3.0.1.sorted.chrom.sizes \
100 \
~/refseqs/dummy_blacklist.bed \
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
samples_rename_merge.txt \
3:00:00 \
20 \
trimmed_align_bowtie2 \
deeptools-hacked_v3.8 \
normal_res \
~/refseqs/dummy_blacklist.bed \
a_agfs_ps


rsync -r /scratch/user/uqgventu/analysis/trimmed_align_* /QRISdata/Q9486/data/UMR_test/

