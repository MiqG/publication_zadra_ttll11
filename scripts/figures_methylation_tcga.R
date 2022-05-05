
require(tidyverse)
require(ggpubr)
require(cowplot)
require(writexl)
require(scattermore)
require(ggrepel)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
ORDER_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP',
             'THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
TEST_METHOD = "wilcox.test"
GENE_OI = "TTLL11"

# formatting
FONT_SIZE = 7 # pt
FONT_FAMILY = 'Helvetica'
PAL_SAMPLE_TYPE = rev(get_palette("npg",2))

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
methylation_file = file.path(PREP_DIR,"methylation_TCGA-tss-by_gene",'merged.tsv.gz')
genexpr_file = file.path(PREP_DIR,'genexpr_TTLLs.tsv')

# outputs
figs_dir = file.path(RESULTS_DIR,"figures","methylation_tcga")

# load
metadata = read_tsv(phenotype_file)
methylation = read_tsv(methylation_file)
genexpr = read_tsv(genexpr_file) %>% filter(gene %in% GENE_OI)

plts = list()

# is TTLL11 differentially methylated across cancers?
X = methylation %>%
    filter(symbol %in% GENE_OI) %>%
    pivot_longer(-c(ensembl_gene,symbol), 
                 names_to="sample", 
                 values_to="beta") %>%
    left_join(metadata, by="sample") %>%
    mutate(sample_type = factor(sample_type, levels=c("Solid Tissue Normal", "Primary Tumor")))
    
plts[["methylation_tcga-TTLL11-counts"]] = X %>%
    count(cancer_type, sample_type) %>%
    ggbarplot(x="cancer_type", y="n", position=position_dodge(0.7),
              fill="sample_type", palette=PAL_SAMPLE_TYPE,
              color=NA, order=ORDER_OI,
              label=TRUE, lab.size=2, lab.family=FONT_FAMILY) +
    labs(x="Cancer Type", y="N. Samples", fill="Sample Type")

plts[["methylation_tcga-TTLL11-boxplots"]] = X %>%
    ggboxplot(x="cancer_type", y="beta", 
              fill="sample_type", palette=PAL_SAMPLE_TYPE,
              order=ORDER_OI, outlier.size = 0.1, lwd=0.1) +
    stat_compare_means(aes(group=sample_type), 
                       method=TEST_METHOD, label="p.signif",
                       size=2, family=FONT_FAMILY) +
    labs(x="Cancer Type", y="median(Methylation Beta)", 
         fill="Sample Type")

methylation_TTLL11 = X

# does the expression of TTLL11 correlate with its methylation?
X = genexpr %>%
    left_join(
        methylation_TTLL11, 
        by=c("gene"="symbol","sample","cancer_type","sample_type")
    ) %>%
    mutate(sample_type = factor(sample_type, levels=c("Solid Tissue Normal", "Primary Tumor")))

plts[["methylation_tcga-cor_expression_met-TTLL11"]] = X %>% 
    ggplot(aes(x=beta, y=expression)) + 
    geom_scattermore(aes(color=sample_type), pixels=c(1000,1000), 
                     pointsize=4, alpha=0.5) + 
    facet_wrap(~sample_type, scales="free_y") + 
    color_palette(PAL_SAMPLE_TYPE) + 
    theme_pubr() + 
    theme() + 
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    labs(x="median(Methylation Beta)", y="Norm. Gene Expression",
         color="Sample Type") +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY),
          aspect.ratio=1)

# does the expresson of TTLL11 correlate with methylation of any other gene?
tumor_samples = genexpr %>% filter(sample_type=="Primary Tumor") %>% pull(sample) %>% unique()
common_samples = intersect(colnames(methylation), tumor_samples)

genexpr_oi = genexpr %>% dplyr::select(sample,expression) %>% deframe()
corr = cor(methylation[,common_samples] %>% t(), 
           genexpr_oi[common_samples], 
           method="spearman", use="pairwise.complete.obs") %>%
    as.vector()

corr = data.frame(
    symbol = methylation[["symbol"]],
    correlation = corr
)

X = corr %>% mutate(gene = "TTLL11")

plts[["methylation_tcga-cor_expression_met-TTLL11_vs_all"]] = X %>% 
    ggviolin(x="gene", y="correlation", color=NA, fill="orange") + 
    geom_boxplot(width=0.1, outlier.size=0.1) + 
    labs(x="Gene", y="Correlation") + 
    geom_text_repel(aes(label=symbol), 
                    X %>% slice_max(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY) + 
    geom_text_repel(aes(label=symbol), 
                    X %>% slice_min(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY)

# save
save_plt = function(plts, plt_name, extension='.pdf', 
                    directory='', dpi=350, format=TRUE,
                    width = par("din")[1], height = par("din")[2]){
    print(plt_name)
    plt = plts[[plt_name]]
    if (format){
        plt = ggpar(plt, font.title=8, font.subtitle=8, font.caption=8, 
                    font.x=8, font.y=8, font.legend=6,
                    font.tickslab=6, font.family=FONT_FAMILY)   
    }
    filename = file.path(directory,paste0(plt_name,extension))
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm', device=cairo_pdf)
}

## plots
dir.create(figs_dir, recursive=TRUE)
save_plt(plts, "methylation_tcga-TTLL11-counts", ".pdf", figs_dir, width=12, height=7)
save_plt(plts, "methylation_tcga-TTLL11-boxplots", ".pdf", figs_dir, width=12, height=7)
save_plt(plts, "methylation_tcga-cor_expression_met-TTLL11", ".pdf", figs_dir, width=8, height=8)
save_plt(plts, "methylation_tcga-cor_expression_met-TTLL11_vs_all", ".pdf", figs_dir, width=6, height=6)

## figdata
