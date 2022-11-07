# Chromosome segregation fidelity requires microtubule polyglutamylation by the cancer downregulated enzyme TTLL11

This repository contains all the scripts to perform the bioinformatic data analyses in Zadra *et al.*, XXXX (DOI: XXXX).

## Overview
The scripts carry out the following steps from data download to publishable figures:

0. **Download data**
1. **Preprocess data**
2. **Expression of TTLL11 and TTLL13 across human tissues.**
3. **Differential expression of TTLL11 across cancer types.**
4. **Correlation between aneuploidy score and expression of TTLL11 across cancer types.**
5. **Mutation frequency per kilobase of every gene in primary tumors.**
6. **Co-expression of TTLLs in healthy and tumor tissues**
7. **Differential expression in cell lines over-expression oncogenes**
8. **Promoter methylation of TTLL11 in tumors**
9. **Differential expression of TTLL11 and TTLL13 in tumors**

## Requirements
In brackets, the version used in this project.
- snakemake (5.31.1)
- gtftools (0.8.5)
- Python (3.8.3)
    - pandas (1.3.0)
- R (4.1.0)
    - clusterProfiler (4.0.5)
    - ComplexHeatmap (2.9.3)
    - cowplot (1.1.1)
    - doParallel (1.0.16)
    - enrichplot (1.12.2)
    - extrafont (0.17)
    - GeneBreak (1.22.0)
    - GenomicRanges (1.44.0)
    - GEOquery (2.60.0)
    - ggplotify (0.1.0)
    - ggpubr (0.4.0)
    - ggrepel (0.9.1)
    - ggvenn (0.1.9)
    - gtools (3.9.2)
    - latex2exp (0.5.0)
    - limma (3.48.3)
    - magrittr (2.0.3)
    - optparse (1.7.1)
    - org.Hs.eg.db (3.13.0)
    - patchwork (1.1.1)
    - pheatmap (1.0.12)
    - readxl (1.3.1)
    - reshape2 (1.4.4)
    - scattermore (0.7)
    - tidyverse (1.3.1)
    - writexl (1.4.0)

## Usage
```shell
git clone https://github.com/MiqG/publication_zadra_ttll11
cd publication_zadra_ttll11
snakemake --cores 10
```

