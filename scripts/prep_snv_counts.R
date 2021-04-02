# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Count mutation types per gene per sample.
#

require(readr)
require(dplyr)
require(tidyr)
require(GeneBreak)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')
RESULTS_DIR = file.path(ROOT,'results')

# inputs
raw_snv_file = file.path(TCGA_DIR,'snv','mc3.v0.2.8.PUBLIC.xena.gz')
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')

# outputs
snv_counts_file = file.path(RESULTS_DIR,'files','snv_gene_freq.tsv')

# load data
## phenotype
phenotype = read_tsv(phenotype_file)
## mutations
snv = read_tsv(raw_snv_file)
## genome annotations
data( ens.gene.ann.hg38 )
genome = ens.gene.ann.hg38
genome = genome %>% mutate(length = End - Start)

# add cancer type
snv = merge(snv, phenotype, by='sample')

# count mutations per gene per cancer type
snv = snv %>% count(cancer_type, gene, effect)

# correct to frequency per kilobase
snv$gene_length = genome$length[match(snv$gene,genome$Gene)]
snv$freq_per_kb = (snv$n / snv$gene_length) * 1000

# remove NAs
snv = snv %>% drop_na()

# save
write_tsv(snv, snv_counts_file)


print('Done!')