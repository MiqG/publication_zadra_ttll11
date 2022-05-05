#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Subset expression of TTLL11 across samples and combine
# with sample metadata.
#

require(tidyverse)
require(limma)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')

# variables
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(TCGA_DIR,'rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')

# outputs
output_file = file.path(PREP_DIR,'genexpr_TCGA.tsv.gz')

# load data
metadata = read_tsv(phenotype_file) %>% 
            filter(cancer_type %in% CANCERS_OI) %>%
            filter(sample_type %in% SAMPLE_TYPES_OI)
samples_oi = metadata %>% pull(sample)
mat = read.columns(genexpr_file, required.col = c('sample',samples_oi))

# preprocess
genexpr = mat 

# save
write_tsv(genexpr, output_file)

print('Done!')
