#!/bin/bash
# Retrive, merge, and relabel V13 forward and reverse reads based in a sample list (samples.txt) and merge them into one file.

mkdir V13_raw_data

while read SAMPLES
do
NAME=$SAMPLES;
find Sequence_data/ -name $NAME\_*R1* -exec gzip -cd {} \; > V13_raw_data/$NAME\_R1.fq 
find Sequence_data/ -name $NAME\_*R2* -exec gzip -cd {} \; > V13_raw_data/$NAME\_R2.fq 
usearch11 -fastq_mergepairs V13_raw_data/$NAME\_R1.fq -sample @ -fastqout V13_raw_data/temp.$NAME.fq
usearch11 -fastx_relabel V13_raw_data/temp.$NAME.fq -prefix $NAME\_ -fastqout V13_raw_data/temp2.$NAME.fq
cat V13_raw_data/temp2.$NAME.fq >> V13_raw_data/V13_merged.fastq
done < V13_samples.txt
date

rm V13_raw_data/temp*.fq

# Usearch pipeline

mkdir V13_Usearch_analysis

### Quality filtering
usearch11 -filter_phix V13_raw_data/V13_merged.fastq -output V13_Usearch_analysis/V13_phixfiltered.fq
usearch11 -fastq_filter V13_Usearch_analysis/V13_phixfiltered.fq -fastq_maxee 1.0 -fastaout V13_Usearch_analysis/V13_qcfiltered.fa

### Identify unique sequences
usearch11 -fastx_uniques V13_Usearch_analysis/V13_qcfiltered.fa -fastaout V13_Usearch_analysis/V13_unique.fa -sizeout -relabel uniq

### Create ASVs
usearch11 -unoise3 V13_Usearch_analysis/V13_unique.fa -zotus V13_Usearch_analysis/V13_ASV.fa
sed -i 's/Zotu/ASV/g' V13_Usearch_analysis/V13_ASV.fa 

### Create ASV table
usearch11 -otutab V13_Usearch_analysis/V13_qcfiltered.fa -zotus V13_Usearch_analysis/V13_ASV.fa -sample_delim _\
     -strand plus -otutabout V13_Usearch_analysis/V13_ASVtab.txt -threads 10

### Predict taxonomy using sintax (MiDAS5.1)
usearch11 -sintax V13_Usearch_analysis/V13_ASV.fa -db /user_data/md/Projects/P026/MiDAS5/MiDAS5.1_20230726/output/FLASVs_w_sintax.fa\
     -tabbedout V13_Usearch_analysis/V13_ASV_midas5.1.sintax -strand plus -sintax_cutoff 0.8 -threads 10
