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
require(magrittr)
require(reshape2)
require(ggpubr)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')

# variables
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')
GENES_OI = c('TTLL1', 'TTLL2', 'TTLL4', 'TTLL5', 'TTLL6', 'TTLL7', 'TTLL9', 'TTLL11', 'TTLL13')

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(TCGA_DIR,'rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')

# outputs
output_file = file.path(PREP_DIR,'genexpr_TTLLs.tsv')

# load data
metadata = read_tsv(phenotype_file) %>% 
            filter(cancer_type %in% CANCERS_OI) %>%
            filter(sample_type %in% SAMPLE_TYPES_OI)
samples_oi = metadata %>% pull(sample)
mat = read.columns(genexpr_file, required.col = c('sample',samples_oi))

# preprocess
idx = mat$sample %in% GENES_OI
genexpr = mat[idx,] %>%
            data.frame(row.names=c()) %>%
            column_to_rownames('sample') %>%
            t() %>%
            melt(varnames = c('sample', 'gene'), value.name = 'expression') %>%
            mutate(sample = gsub('\\.','-',sample),
                   expression = as.numeric(expression)) %>%
            left_join(metadata, by='sample') %>%
            drop_na()

# save
write_tsv(genexpr, output_file)

print('Done!')