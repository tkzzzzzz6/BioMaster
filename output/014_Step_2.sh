#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/014
conda install -y -c bioconda star samtools
mkdir -p ./output/014/star_index
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./output/014/star_index --genomeFastaFiles ./data/minigenome.fa --genomeSAindexNbases 8
STAR --runThreadN 4 --genomeDir ./output/014/star_index --readFilesIn ./output/014/trimmed_rnaseq_1.fastq.gz ./output/014/trimmed_rnaseq_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./output/014/ --outSAMtype BAM SortedByCoordinate --outBAMcompression 5
mv ./output/014/Aligned.sortedByCoord.out.bam ./output/014/aligned_rnaseq.bam
samtools index ./output/014/aligned_rnaseq.bam
