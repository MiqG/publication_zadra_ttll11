# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Visualize mutation frequencies across all genes.
#
# Outline
# -------
# - normalized mutation frequencies by cancer and mutation effect
# - median normalized mutation frequencies across cancers by mutation effect


require(readr)
require(dplyr)
require(tidyr)
require(ggpubr)
require(ggrepel)
require(latex2exp)
require(gtools)
require(writexl)

# variables
ROOT = here::here()
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results')

GENE_OI = 'TTLL11'
ORDER_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')

WIDTH = 7.5
HEIGHT = 4
BASE_SIZE = 10

# inputs
snv_counts_file = file.path(PREP_DIR,'snv_gene_freq.tsv')

# outputs
output_file = file.path(RESULTS_DIR,'figures','mutation_frequency.rds')
output_figdata = file.path(RESULTS_DIR,'files','figdata-mutation_frequency.xlsx')


##### FUCNTIONS #####
make_plot_bycancer = function(df, gene_oi=GENE_OI, y_lab_pos=-5, is_pancan=FALSE){
    if(!is_pancan){
        ## add missing cancer in gene OI
        missing_cancers = setdiff(ORDER_OI, df %>% filter(gene == gene_oi) %>% distinct(cancer_type,gene) %>% pull(cancer_type))
        df = df %>% ungroup() %>% add_row(cancer_type = missing_cancers, gene='TTLL11', value=NULL, pvalue=1)

        ## reorder cancer types
        df = df %>% mutate(cancer_type = factor(cancer_type, levels=ORDER_OI))
        df = df[order(df$cancer_type),]


    }
    ## prepare pvalues in gene OI
    stat.test = df %>% filter(gene == gene_oi) %>% 
    mutate(
        p=pvalue, 
        p.format=format.pval(pvalue,1),
        p.signif=stars.pval(pvalue),
        group1=NA,
        group2=NA)

    ## plot
    plt = ggviolin(df, x='cancer_type', y='z_score', color='darkgrey', alpha=0.5, lwd=0.2) +
    geom_boxplot(width=0.1, outlier.shape = NA, lwd=0.1) + 
    geom_point(data = df %>% filter(gene == 'TTLL11'), color='darkred', size=1) + 
    xlab(element_blank()) + ylab(element_blank()) + ggtitle(element_blank()) + 
    stat_pvalue_manual(stat.test, x = 'cancer_type', y.position = y_lab_pos, label = 'p.signif')
    plt = ggpar(plt, font.tickslab = c(BASE_SIZE))
    return(plt)
}


##### SCRIPT #####
# load data
snv = read_tsv(snv_counts_file)

# prepare
## keep only cancer types of interest
snv = snv %>% filter(cancer_type %in% ORDER_OI)

## keep only those effects and cancers where TTLL11 appears
effects_oi = snv %>% filter(gene %in% GENE_OI) %>% distinct(effect) %>% pull()
snv = snv %>% filter(effect %in% effects_oi)


##### Plot mutation frequency by cancer and effect #####
df = snv

# compute z-score and pvalue by cancer and effect
df = df %>% 
    group_by(cancer_type,effect) %>% 
    mutate(
        z_score=as.numeric(scale(log10(freq_per_kb))),
        pvalue=pnorm(-abs(z_score))
    )

# plot
plts = list()
for (eff in unique(df$effect)){
    plt_name = paste0('effect_by_cancer-',eff)
    
    # by cancer
    plt_bycancer = make_plot_bycancer(df %>% filter(effect==eff))
    
    # pancancer
    df_pancancer = df %>% 
        group_by(gene,effect) %>% 
        summarize(freq_per_kb=median(freq_per_kb))
    df_pancancer = df_pancancer %>% 
        group_by(effect) %>% 
        mutate(z_score=as.numeric(scale(log10(freq_per_kb))),
               pvalue=pnorm(-abs(z_score)))
    df_pancancer$cancer_type = 'PANCAN'
    plt_pancancer = make_plot_bycancer(df_pancancer %>% filter(effect==eff), is_pancan = TRUE)
    
    # combine
    #plt = ggarrange(plt_bycancer, plt_pancancer, ncol=2, widths=c(1.5,0.3))
    #plt = annotate_figure(plt,
    #                      top = paste('Gene Mutation Frequency per Kb -',eff),
    #                      left = text_grob(TeX('Normalized $log_{10}$(Mut. Freq.)'), rot = 90),
    #                      bottom = 'Cancer Type')
    
    plts[[plt_name]] = list(bycancer=plt_bycancer, pancancer=plt_pancancer)
}

# prepare figure data
df_pancan = df %>% 
    group_by(gene, gene_length, effect) %>% 
    summarize(n=n(), freq_per_kb=median(freq_per_kb)) %>% 
    ungroup() %>% 
    mutate(z_score=as.numeric(scale(freq_per_kb)), 
           pvalue=pnorm(-abs(z_score))) %>% 
    mutate(cancer_type='PANCAN')

common_cols = intersect(colnames(df),colnames(df_pancan))
figdata = rbind(df[,common_cols], df_pancan[,common_cols])

# save plots
saveRDS(plts, output_file)
## figure data
write_xlsx(
    x = list(mutation_frequency = figdata),
    path = output_figdata
)

print('Done!')