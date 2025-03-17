#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -y -c bioconda plink
mkdir -p ./output/002
plink --bfile ./data/1000GP_pruned --homozyg --homozyg-snp 50 --homozyg-kb 1000 --homozyg-density 50 --homozyg-gap 100 --out ./output/002/roh_results
