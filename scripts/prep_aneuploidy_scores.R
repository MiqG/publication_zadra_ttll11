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
require(readxl)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
ARTICLES_DIR = file.path(DATA_DIR,'raw','articles')
PREP_DIR = file.path(DATA_DIR,'prep')

# inputs
raw_aneuploidy_scores_file = file.path(ARTICLES_DIR,'Taylor2018','aneuploidy.xlsx')

# outputs
prep_aneuploidy_scores_file = file.path(PREP_DIR,'aneuploidy.tsv')

# load data
df = read_excel(raw_aneuploidy_scores_file, sheet = 1, skip = 1)

# subset
df = df[,c('Sample','AneuploidyScore(AS)')]

# rename
colnames(df) = c('sample','aneuploidy_score')

# save
write_tsv(df, prep_aneuploidy_scores_file)


print('Done!')