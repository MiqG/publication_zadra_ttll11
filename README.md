# Microtubule polyglutamylation by TTLL11 ensures mitotic chromosome segregation fidelity

This repository contains all the scripts to perform the bioinformatic data analyses in Zadra *et al.* XXXX (DOI: XXXX).

## Overview
The scripts carry out the following steps:
0. Download data
1. Preprocess data:
2. Expression of TTLL11 and TTLL13 across human tissues.
3. Differential expression of TTLL11 across cancer types.
4. Correlation between aneuploidy score and expression of TTLL11 across cancer types.
5. Mutation frequency per kilobase of every gene in primary tumors.
6. Prepare publishable figures

## Analysis
- [ ] coexpression enrichment of transcription factors [CHEA harmonizome](https://maayanlab.cloud/Harmonizome/dataset/CHEA+Transcription+Factor+Targets) for all TTLLs
- [ ] Do we see the same causality of the overexpressed genes with clear direction in cancers?

## Requirements
In brackets, the version used in this project.
- python 3 (3.8.3)
- snakemake (5.31.1)
- R (4.0.3)
    - tidyverse
    - ggpubr
    - pheatmap
    - doParallel
    - writexl
    - reshape2
    - magrittr
    - latex2exp
    - GeneBreak
    - patchwork
    - gtools

## Usage
```shell
git clone https://github.com/MiqG/Zadra_2021.git
cd Zadra_2021
snakemake --cores 10
```

