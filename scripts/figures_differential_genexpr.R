# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compare gene expression of TTLL11 in TCGA.
#

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
genexpr_file = file.path(RESULTS_DIR,'files','genexpr_TTLL11.tsv')

# outputs
output_file = file.path(RESULTS_DIR,'figures','differential_expression.rds')

# load data
df = read_tsv(genexpr_file)

# change sample_type and cancer_type order
df = df %>% 
    mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPES_OI),
           cancer_type = factor(cancer_type, levels = CANCERS_OI))

# plot
make_plot = function(X){
    palette = rev(get_palette('npg',2))
    plt = ggplot(X, aes(cancer_type, expression, fill=sample_type)) + 
        geom_boxplot(outlier.size = 0.1, lwd=0.1) +
        theme_pubr() + 
        stat_compare_means(aes(group = sample_type), 
                           method=TEST_METHOD, label = "p.signif")
    
    plt = set_palette(plt, palette = palette) +
        labs(x=element_blank(), y=element_blank(), fill='Sample Type')
    return(plt)
}

plts = list()
plts[['bycancer']] = make_plot(df) + scale_y_continuous(limits = quantile(df$expression, c(0.001, 0.999))) 
plts[['pancancer']] = make_plot(df %>% mutate(cancer_type=factor('PANCAN'))) + scale_y_continuous(limits = quantile(df$expression, c(0.001, 0.999)))

# save
saveRDS(plts, output_file)


print('Done!')