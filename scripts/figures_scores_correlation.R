# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Make figure of correlations between aneuploidy and centromere amplification scores 
# with gene expression by cancer and sample types.
#

# libraries
require(readr)
require(ggpubr)
require(dplyr)
require(tidyr)
require(reshape2)
require(latex2exp)
require(gtools)

GENE_OI = 'TTLL11'
ORDER_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')

# variables
ROOT = here::here()
RESULTS_DIR = file.path(ROOT,'results')

WIDTH = 7.5
HEIGHT = 4
BASE_SIZE = 10

N_SAMPLE = 1000

# inputs
correlations_ca_file = file.path(RESULTS_DIR,'files','correlation-genexpr_centrosome_amplification.tsv')
correlations_aneuploidy_file = file.path(RESULTS_DIR,'files','correlation-genexpr_aneuploidy.tsv') 

# outputs
output_file = file.path(RESULTS_DIR,'figures','correlation_with_scores.rds')


##### FUNCTIONS ######
.make_plot = function(df, plt_title, gene_oi=GENE_OI, y_lab_pos=-0.6){
    # plot spearman correlation values as stripchart and violin
    # highlighting TTLL11 and adding the p-value of the Z-score of 
    # TTLL11 as labels
    
    ## prepare pvalues in gene OI
    stat.test = df %>% filter(gene == GENE_OI) %>% 
        mutate(
            p=pvalue, 
            p.format=format.pval(pvalue,1),
            p.signif=stars.pval(pvalue),
            group1=NA,
            group2=NA)

    ## plot
    plt = ggstripchart(df, x='cancer_type', y='value', color='darkgrey', alpha=.1, size=0.5) +
    geom_violin(alpha=0.5, lwd=0.2) + 
    geom_boxplot(width=0.1, outlier.shape = NA, lwd=0.1) + 
    geom_point(data = df %>% filter(gene == 'TTLL11'), color='darkred', size=1) + 
    xlab(element_blank()) + ylab(element_blank()) + ggtitle(element_blank()) + 
    stat_pvalue_manual(stat.test, x = 'cancer_type', y.position = y_lab_pos, label = 'p.signif')
    plt = ggpar(plt, font.tickslab = c(BASE_SIZE))
    return(plt)
}

make_plot = function(df, plt_title, gene_oi=GENE_OI, y_lab_pos=-0.6){
    # plot spearman correlation values as stripchart and violin
    # highlighting TTLL11 and adding the p-value of the Z-score of 
    # TTLL11 as labels
    
    ## prepare pvalues in gene OI
    stat.test = df %>% filter(gene == GENE_OI) %>% 
        mutate(
            p=pvalue, 
            p.format=format.pval(pvalue,1),
            p.signif=stars.pval(pvalue),
            group1=NA,
            group2=NA)

    ## plot
    plt = ggviolin(df, x='cancer_type', y='value', color='darkgrey', alpha=0.5, lwd=0.2) +
    geom_boxplot(width=0.1, outlier.shape = NA, lwd=0.1) + 
    geom_point(data = df %>% filter(gene == 'TTLL11'), color='darkred', size=1) + 
    xlab(element_blank()) + ylab(element_blank()) + ggtitle(element_blank()) + 
    stat_pvalue_manual(stat.test, x = 'cancer_type', y.position = y_lab_pos, label = 'p.signif')
    plt = ggpar(plt, font.tickslab = c(BASE_SIZE))
    return(plt)
}

##### Data wrangling ######
# load data
correlations_ca = read_tsv(correlations_ca_file)
correlations_aneuploidy = read_tsv(correlations_aneuploidy_file)

# transform to long format
correlations_ca = correlations_ca %>% melt() %>% separate(variable,c('cancer_type','sample_type')) %>% mutate(score_type = 'Centrosome Amplification (CA20)', sample_type = gsub("([a-z])([A-Z])","\\1 \\2", sample_type))
correlations_aneuploidy = correlations_aneuploidy %>% melt() %>% rename(cancer_type = variable) %>% mutate(sample_type = 'Primary Tumor', score_type = 'Aneuploidy') 


# filter
correlations_ca = correlations_ca %>% filter(cancer_type %in% ORDER_OI)
correlations_aneuploidy = correlations_aneuploidy %>% filter(cancer_type %in% ORDER_OI)

plts = list()

##### PLOTS ANEUPLOIDY #####
# init
df = correlations_aneuploidy

# compute z-scores and pvalues by cancer
df = df %>% group_by(cancer_type) %>% mutate(z_score=scale(value),pvalue=pnorm(-abs(z_score))) %>% ungroup()

# change group order
df = df %>% mutate(cancer_type = factor(cancer_type, levels=ORDER_OI))
df = df[order(df$cancer_type),]

# plot by cancer
df = df %>% drop_na()
plt_title = 'Gene Expression and Aneuploidy Score'
gene_oi = GENE_OI
plts[['aneuploidy_bycancer']] = make_plot(df, plt_title, gene_oi, -0.7)

# prepare pancancer
df_pancan = df %>% group_by(gene) %>% summarize(value=median(value)) %>% ungroup() %>% mutate(z_score=as.numeric(scale(value)), pvalue=pnorm(-abs(z_score))) %>% mutate(cancer_type='PANCAN')

# plot pancancer
plts[['aneuploidy_pancancer']] = make_plot(df_pancan, plt_title, gene_oi, -0.4)

# save
saveRDS(plts, output_file)

print('Done!')