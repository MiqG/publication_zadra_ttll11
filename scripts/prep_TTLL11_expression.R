require(tidyverse)
require(limma)
require(magrittr)
require(ggpubr)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')
GENE_OI = 'TTLL11'

TEST_METHOD = 'wilcox.test'

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(TCGA_DIR,'rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')

# outputs
output_file = file.path(RESULTS_DIR,'files','genexpr_TTLL11.tsv')

# load data
metadata = read_tsv(phenotype_file) %>% 
            filter(cancer_type %in% CANCERS_OI) %>%
            filter(sample_type %in% SAMPLE_TYPES_OI)
samples_oi = metadata %>% pull(sample)
mat = read.columns(genexpr_file, required.col = c('sample',samples_oi))

# preprocess
idx = mat$sample %in% GENE_OI
genexpr = t(mat[idx,]) %>% 
            as.data.frame() %>% 
            rownames_to_column() %>% 
            set_colnames(c('sample', 'expression')) %>% 
            filter(expression!=GENE_OI) %>% 
            mutate(sample = gsub('\\.','-',sample),
                   expression = as.numeric(expression))

# combine
df = merge(metadata, genexpr, by='sample')

# save
write_tsv(df, output_file)

print('Done!')