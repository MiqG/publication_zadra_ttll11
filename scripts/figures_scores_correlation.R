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

GENE_OI = 'TTLL11'
ORDER_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')

# variables
ROOT = here::here()
RESULTS_DIR = file.path(ROOT,'results')

WIDTH = 6.5
HEIGHT = 4
BASE_SIZE = 10

# inputs
correlations_ca_file = file.path(RESULTS_DIR,'files','correlation-genexpr_centrosome_amplification.tsv')
correlations_aneuploidy_file = file.path(RESULTS_DIR,'files','correlation-genexpr_aneuploidy.tsv') 

# outputs
fig_dir = file.path(RESULTS_DIR,'figures','correlation_with_scores')

##### FUNCTIONS ######
sorted_stripchart = function(df, plt_title, genes_oi=GENE_OI, order_oi=ORDER_OI){
    ## change group order
    df = df %>% mutate(cancer_type = factor(cancer_type, levels=order_oi))
    df = df[order(df$cancer_type),]
    ## sort and add id
    df = df %>% drop_na() %>% arrange(cancer_type,sample_type,value)
    df$id = 1:nrow(df)
    ## make scatter
    plt = df %>% ggplot(aes(x=id, y=value)) + geom_point(alpha=0.5, color='darkgray') + geom_point(data=df %>% filter(gene %in% genes_oi), aes(x=id, y=value), color='darkred', size=2) + geom_hline(yintercept = c(-1.96,1.96), linetype='dashed')
    ## change ticks
    tick_breaks = df %>% group_by(cancer_type) %>% mutate(tick_breaks=median(id)) %>% distinct(tick_breaks) %>% pull()
    tick_labels = df %>% distinct(cancer_type) %>% pull()
    plt = plt + scale_x_continuous(breaks = tick_breaks, labels = tick_labels) 
    ## make pretty
    plt = plt + theme_pubr(base_size=BASE_SIZE)
    plt = plt + ylab('Normalized Spearman Correlation (Z-score)') + xlab('Cancer Type') + ggtitle(plt_title)
    
    return(plt)
}

# load data
correlations_ca = read_tsv(correlations_ca_file)
correlations_aneuploidy = read_tsv(correlations_aneuploidy_file)

# transform to long format
correlations_ca = correlations_ca %>% melt() %>% separate(variable,c('cancer_type','sample_type')) %>% mutate(score_type = 'Centrosome Amplification (CA20)', sample_type = gsub("([a-z])([A-Z])","\\1 \\2", sample_type))
correlations_aneuploidy = correlations_aneuploidy %>% melt() %>% rename(cancer_type = variable) %>% mutate(sample_type = 'Primary Tumor', score_type = 'Aneuploidy') 

# filter
correlations_ca = correlations_ca %>% filter(cancer_type %in% ORDER_OI)
correlations_aneuploidy = correlations_aneuploidy %>% filter(cancer_type %in% ORDER_OI)

# plot overall distributions
plts = list()
plts[['aneuploidy']] = sorted_stripchart(correlations_aneuploidy %>% group_by(cancer_type) %>% mutate(value=scale(value)), 'Gene Expression and Aneuploidy Score')
plts[['centrosome_amplification_pt']] = sorted_stripchart(correlations_ca %>% filter(sample_type=='Primary Tumor') %>% mutate(value=scale(value)), 'Gene Expression and Centrosome Amplification Score in PT')
plts[['centrosome_amplification_stn']] = sorted_stripchart(correlations_ca %>% filter(sample_type=='Solid Tissue Normal') %>% mutate(value=scale(value)), 'Gene Expression and Centrosome Amplification Score in STN')

# save
lapply(names(plts), function(plt_name){ 
    filename = file.path(fig_dir,paste0(plt_name,'.png'))
    ggsave(filename, plts[[plt_name]], width = WIDTH, height = HEIGHT, dpi = 100) 
})

# summary of TTLL11 across cancers
df = rbind(correlations_ca, correlations_aneuploidy[,colnames(correlations_ca)])
plt = df %>% filter(gene == GENE_OI) %>% ggstripchart(x = 'sample_type', y = 'value', color='cancer_type', palette = get_palette('Spectral',length(ORDER_OI))) + facet_wrap(~score_type, scales='free') + ylab('Spearman Correlation') + xlab('Sample Type') + ggtitle(TeX(r'(\textit{TTLL11})')) + labs(color = 'Cancer Type') + theme_pubr(base_size=BASE_SIZE, legend = 'bottom')

#df %>% filter(gene == GENE_OI & score_type=='Centrosome Amplification (CA20)') %>% ggpaired(x = 'sample_type', y = 'value', color='sample_type', line.color = "gray", line.size = 0.4) 

# save
filename = file.path(fig_dir,'correlations_TTLL11.pdf')
ggsave(filename, plt, dpi=100, height = HEIGHT, width = WIDTH)

print('Done!')