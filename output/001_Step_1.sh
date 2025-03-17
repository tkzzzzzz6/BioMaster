#!/bin/bash
which python
conda config --set show_channel_urls false
conda config --add channels conda-forge
conda config --add channels bioconda
mkdir -p ./output/001
conda install -y plink
plink --bfile ./data/1000GP_pruned --r2 --ld-window 99999 --ld-window-kb 500 --out ./output/001/ld_decay
