#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/005
conda install -y trimmomatic
trimmomatic PE -phred33 ./data/rnaseq_1.fastq.gz ./data/rnaseq_2.fastq.gz ./output/005/trimmed_rnaseq_1.fastq.gz ./output/005/unpaired_rnaseq_1.fastq.gz ./output/005/trimmed_rnaseq_2.fastq.gz ./output/005/unpaired_rnaseq_2.fastq.gz ILLUMINACLIP:./data/minigenome.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
