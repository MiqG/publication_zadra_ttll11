# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Standardize column names and clean table from article. 
#

# libraries
require(readr)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
ARTICLES_DIR = file.path(DATA_DIR,'raw','articles')
PREP_DIR = file.path(DATA_DIR,'prep')

# inputs
raw_ca_scores_file = file.path(ARTICLES_DIR,'Almeida2019','centrosome_amplification.tsv')

# outputs
prep_ca_scores_file = file.path(PREP_DIR,'centrosome_amplification.tsv')

# load data
df = read_tsv(raw_ca_scores_file)

# subset
df = df[,c('Sample ID','CA20')]

# rename
colnames(df) = c('sample','CA20')

# reformat sample barcodes
df$sample = gsub('\\.','-',df$sample)

# save
write_tsv(df, prep_ca_scores_file)


print('Done!')