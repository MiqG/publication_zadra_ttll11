# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Standardize column names and clean table. 
# Keep only information for those cancer types that have >20 STN and PT samples.
#

# libraries
require(readr)
require(dplyr)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')

SAMPLE_TYPES_OI = c('Primary Tumor','Solid Tissue Normal')
THRESH = 20

# inputs
raw_survival_file = file.path(TCGA_DIR,'phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz')
raw_phenotype_file = file.path(TCGA_DIR,'phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz')

# outputs
prep_phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')

# load data
df = read_tsv(raw_phenotype_file)
df_surv = read_tsv(raw_survival_file)

# add cancer type abbreviations to phenotype
abbreviations = merge(df, df_surv, by='sample') %>% distinct(`cancer type abbreviation`,`_primary_disease`)
df = merge(df, abbreviations, by='_primary_disease')

# subset
df = df[,c('sample','cancer type abbreviation','sample_type')]

# rename
colnames(df) = c('sample','cancer_type','sample_type')

# select samples and cancer types of interest
## which cancer types have more than 'THRESH' samples of every sample type?
cancer_types_oi = df %>% 
                    filter( sample_type %in% SAMPLE_TYPES_OI ) %>%
                    group_by(cancer_type, sample_type) %>% 
                    summarize(n = n()) %>%
                    filter( n > THRESH ) %>%
                    filter( duplicated(cancer_type) ) %>%
                    pull(cancer_type)
## keep only these
df = df %>% filter( sample_type %in% SAMPLE_TYPES_OI & cancer_type %in% cancer_types_oi )

# save
write_tsv(df, prep_phenotype_file)

print('Done!')