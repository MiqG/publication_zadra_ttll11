# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compare gene expression of TTLLs in TCGA.
#

require(tidyverse)
require(limma)
require(magrittr)
require(ggpubr)
require(writexl)
require(latex2exp)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP',
               'THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')

TEST_METHOD = 'wilcox.test'

FONT_SIZE = 7 # pt
FONT_FAMILY = 'helvetica'

# inputs
genexpr_file = file.path(PREP_DIR,'genexpr_TTLLs.tsv')

# outputs
diffexpr_TTLL11_file = file.path(RESULTS_DIR,'figures','differential_expression.pdf')
diffexpr_TTLLs_file = file.path(RESULTS_DIR,'figures','TTLLs-differential_gene_expression-barplot.pdf')
output_figdata = file.path(RESULTS_DIR,'files','figdata-differential_expression.xlsx')

# load data
df = read_tsv(genexpr_file)

# change sample_type and cancer_type order
df = df %>% 
    mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPES_OI),
           cancer_type = factor(cancer_type, levels = CANCERS_OI))

##### Run Differential analysis for all TTLLs across cancers #####
# differential gene expression function
run_wilcox = function(df, condition_a='Primary Tumor', condition_b='Solid Tissue Normal'){
    genes_oi = df %>% pull(gene) %>% unique()
    result = lapply(genes_oi, function(gene_oi){
        x = df %>% 
            filter(gene%in%gene_oi & sample_type%in%condition_a) %>% 
            pull(expression) %>% na.omit()
        y = df %>% 
            filter(gene%in%gene_oi & sample_type%in%condition_b) %>% 
            pull(expression) %>% na.omit()
        
        test = wilcox.test(x,y)
        
        result = data.frame(
            gene = gene_oi,
            logFC = median(x) - median(y),
            pvalue = test[['p.value']],
            fdr = p.adjust(test[['p.value']], method='fdr'),
            condition_a = condition_a,
            condition_b = condition_b,
            n_condition_a = length(x),
            n_condition_b = length(y)
        )
        return(result)
    })
    result = do.call(rbind,result)
    
    return(result)
}

# run analysis
result = lapply(CANCERS_OI, function(cancer_type_oi){
        tmp = df %>% filter(cancer_type %in% cancer_type_oi)
        result = run_wilcox(tmp) %>% mutate(cancer_type = cancer_type_oi)
        return(result)
}) 
result = do.call(rbind, result)

# prepare data
sample_summary = df %>%
    group_by(cancer_type, sample_type, gene) %>%
    summarize(n = n(),
              median = median(expression),
              mean = mean(expression),
              std = sd(expression),
              mad = mad(expression))

# save figure data
write_xlsx(
    x = list(sample_summary = sample_summary, test_result = result),
    path = output_figdata
)

##### Plot summary of differential gene expression analysis of TTLLs #####
# plotting function
plot_summary_dge = function(result){
    
    df = result %>% 
            group_by(cancer_type) %>%
            mutate(log10_fdr = -log10(fdr)) %>% 
            ungroup()
    
    X = df %>% 
        mutate(is_sign = fdr < 0.05, 
               regulation = ifelse(logFC>0,'upregulated','downregulated')) %>% 
        filter(is_sign) %>% 
        count(gene,regulation) %>% 
        mutate(nn = ifelse(regulation=='downregulated',-n,n)) 
    
    genes_order = X %>% 
        group_by(gene) %>% 
        summarise(s=sum(n)) %>% 
        arrange(s) %>% 
        pull(gene)

    plt = X %>% 
        ggbarplot(x='gene', y='nn', fill='regulation', 
                  color=NA, order=genes_order, orientation = "horiz", 
                  palette='simpsons') + geom_text(aes(label=n)) + 
        geom_hline(yintercept=0) + 
        labs(x='Gene', y='Count', fill = 'Regulation') + 
        theme(axis.line = element_blank())
    
    return(plt)
}

# make plot
plt = plot_summary_dge(result)

# save
ggsave(diffexpr_TTLLs_file, plt, units = 'cm', width = 12, height = 12)


##### Plot Gene Expression Distributions of TTLL11 across cancers #####
# plotting function
make_boxplots_gene_oi = function(X, gene_oi){
    palette = rev(get_palette('npg',2))
    plt = X %>%
        filter(gene %in% gene_oi) %>%
        ggplot(aes(cancer_type, expression, fill=sample_type)) + 
        geom_boxplot(outlier.size = 0.1, lwd=0.1) +
        theme_pubr() + 
        stat_compare_means(aes(group = sample_type), 
                           method=TEST_METHOD, label = "p.signif")
    
    plt = set_palette(plt, palette = palette) +
        labs(x=element_blank(), y=element_blank(), fill='Sample Type')
    return(plt)
}

# plot
plts = list()
plts[['bycancer']] = make_boxplots_gene_oi(df,'TTLL11') + scale_y_continuous(limits = quantile(df$expression, c(0.001, 0.999))) 
plts[['pancancer']] = make_boxplots_gene_oi(df %>% mutate(cancer_type=factor('PANCAN')), 'TTLL11') + scale_y_continuous(limits = quantile(df$expression, c(0.001, 0.999)))

# compose figure
widths_cancer_type = c(1.15,0.13)

plt = ggarrange(
    plts[['bycancer']] + 
        ylab(TeX('$log_2(Norm. Count + 1)$')) + 
        theme_pubr(base_size = FONT_SIZE) +
        theme(plot.margin = margin(0,0,0,0, "cm")),
    #NULL,
    plts[['pancancer']] + 
        theme_pubr(base_size = FONT_SIZE) +
        theme(plot.margin = margin(0,0,0,0, "cm")),
    widths = widths_cancer_type,
    common.legend = TRUE,
    ncol = 2
) %>% annotate_figure(bottom = text_grob('Cancer Type', 
                                         size = FONT_SIZE, 
                                         family = FONT_FAMILY,
                                         vjust = -0.5)
                     ) +
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))

# save figure
ggsave(plt, filename=diffexpr_TTLL11_file, units = 'cm', width=12, dpi=300, height=8.5, device=cairo_pdf)


print('Done!')## figure data
write_xlsx(
    x = list(sample_summary = sample_summary, test_result = result),
    path = output_figdata
)