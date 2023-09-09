#!/bin/bash
# Retrive and relabel V4 forward reads based on a sample list (samples.txt) and merge them into one file.

mkdir V4_raw_data

while read SAMPLES
do
NAME=$SAMPLES;
find Sequence_data/ -name $NAME\_*R1* -exec gzip -cd {} \; > V4_raw_data/$NAME.R1.fq 
usearch11 -fastx_relabel V4_raw_data/$NAME.R1.fq -prefix $NAME\_ -fastqout V4_raw_data/temp.$NAME.R1.fq
cat V4_raw_data/temp.$NAME.R1.fq >> V4_raw_data/V4_forward.fastq
done < V4_samples.txt
date

rm V4_raw_data/temp*.fq

# Usearch pipeline
mkdir V4_Usearch_analysis

### Adapter trimming using cutadapt as the exact position of the primers are not known
module load cutadapt/2.8-foss-2018a-Python-3.6.4 # calling in-house installed version of cutadapt
cutadapt -g ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC -o V4_Usearch_analysis/V4_trimmed.fq --discard-untrimmed V4_raw_data/V4_forward.fastq -j 10

### Quality filtering
usearch11 -fastx_revcomp V4_Usearch_analysis/V4_trimmed.fq -fastqout V4_Usearch_analysis/V4_trimmed_rc.fq -threads 10 #reverse complement as our V4 reads are read in reverse
usearch11 -fastq_filter V4_Usearch_analysis/V4_trimmed_rc.fq -fastq_maxee 1.0 -fastaout V4_Usearch_analysis/V4_qcfiltered.fa -threads 10

### Identify unique sequences
usearch11 -fastx_uniques V4_Usearch_analysis/V4_qcfiltered.fa -fastaout V4_Usearch_analysis/V4_unique.fa -sizeout -relabel uniq

### Create ASVs
usearch11 -unoise3 V4_Usearch_analysis/V4_unique.fa -zotus V4_Usearch_analysis/V4_ASV.fa
sed -i 's/Zotu/ASV/g' V4_Usearch_analysis/V4_ASV.fa 

### Create ASV table
usearch11 -otutab V4_Usearch_analysis/V4_qcfiltered.fa -zotus V4_Usearch_analysis/V4_ASV.fa -sample_delim _ -strand plus\
     -otutabout V4_Usearch_analysis/V4_ASVtab.txt -threads 10

### Predict taxonomy using sintax (MiDAS5.1)
usearch11 -sintax V4_Usearch_analysis/V4_ASV.fa -db /user_data/md/Projects/P026/MiDAS5/MiDAS5.1_20230726/output/FLASVs_w_sintax.fa\
     -tabbedout V4_Usearch_analysis/V4_ASV_MiDAS_5.1.sintax -strand plus -sintax_cutoff 0.8 -threads 10
