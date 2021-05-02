# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Measure spearman correlation between aneuploidy scores and gene expression by cancer.
#

# libraries
require(readr)
require(dplyr)
require(limma)
require(tibble)
require(doParallel)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')
RESULTS_DIR = file.path(ROOT,'results')

COR_METHOD = 'spearman'

N_CORES = 10

# inputs
aneuploidy_scores_file = file.path(PREP_DIR,'aneuploidy.tsv')
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(TCGA_DIR,'rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')

# outputs
correlations_file = file.path(RESULTS_DIR,'files','correlation-genexpr_aneuploidy.tsv')

# load data (NOTE: gene expression data is loaded on the fly)
score = read_tsv(aneuploidy_scores_file) 
metadata = read_tsv(phenotype_file)

# available samples 
## in gene expression matrix
genexpr_samples = strsplit(readLines(file(genexpr_file,'r'),n=1),split='\t')[[1]]
genexpr_samples = setdiff(genexpr_samples,'sample')
genexpr_genes = read.columns(genexpr_file, required.col = 'sample')[[1]]

## aneuploidy score
score_samples = score$sample

# measure correlation between gene expression and aneuplo
cancer_types = metadata %>% pull(cancer_type) %>% unique()

cl =  makeCluster(N_CORES)
registerDoParallel(cl)
correlations = foreach(cancer_type_oi=cancer_types, 
                       .combine=cbind,
                       .packages=c('dplyr','limma')) %dopar% {
    # get samples for correlation
    samples_oi = metadata %>% filter(cancer_type==cancer_type_oi) %>% pull(sample)
    samples_oi = intersect(samples_oi, intersect(genexpr_samples, score_samples))

    # read gene expression for these samples
    genexpr = read.columns(genexpr_file, required.col = c('sample',samples_oi))
    genexpr[,samples_oi] = genexpr[,samples_oi]

    # compute correlations
    x = t(genexpr[,samples_oi])
    y = score$aneuploidy_score[match(samples_oi, score$sample)]
    corr = cor(x, y, method=COR_METHOD)    
    colnames(corr) = cancer_type_oi
    return(corr)
}
stopCluster(cl)

# add genes
correlations = data.frame(gene=genexpr_genes, correlations)

# save
write_tsv(correlations, correlations_file)


print('Done!')