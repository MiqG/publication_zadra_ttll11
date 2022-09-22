# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compare TTLL11 and TTLL13
#
# Outline
# -------
# - correlation in tumors
# - distribution of expression of TTLL13 vs TTLL11 in tumors

require(tidyverse)
require(cowplot)
require(ggpubr)
require(writexl)
require(scattermore)
require(extrafont)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
GENES_OI = c("TTLL11","TTLL13")
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP',
               'THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')
TEST_METHOD = "wilcox.test"

# formatting
PAL_GENE = get_palette("Dark2",2)
FONT_SIZE = 7 # pt
FONT_FAMILY = 'Arial'

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(PREP_DIR,'genexpr_TCGA.tsv.gz')

# outputs
figs_dir = file.path(RESULTS_DIR,"figures","TTLL11_vs_TTLL13")

# load
metadata = read_tsv(phenotype_file)
genexpr = read_tsv(genexpr_file)

# prep
metadata = metadata %>%
    filter(cancer_type %in% CANCERS_OI)

genexpr_oi = genexpr %>% 
    dplyr::rename(gene=sample) %>% 
    filter(gene %in% GENES_OI) %>% 
    column_to_rownames("gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample")

X = metadata %>%
    left_join(
        genexpr_oi,
        by = "sample"
    ) %>%
    filter(sample_type=="Primary Tumor") %>%
    mutate(cancer_type = factor(cancer_type, levels=CANCERS_OI))

plts = list()

# does the expresson of TTLL11 correlate with expression of TTLL13?
plts[["coexpression-scatter"]] = X %>%
    ggplot(aes(x=TTLL11, y=TTLL13)) +
    geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
    facet_wrap(~cancer_type) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    theme_pubr() +
    labs(x="NES Primary Tumor", y="NES Solid Tissue Normal")


# how does the expression of TTLL11 differ from TTLL13?
x = X %>%
    pivot_longer(c(TTLL11,TTLL13), names_to="gene", values_to="expression")

tests = compare_means(
    formula=expression~gene, data=x, 
    method=TEST_METHOD, group.by="cancer_type") %>%
    mutate(gene=group2)

plts[["expression-boxplot-bycancer"]] = x %>%
    ggplot(aes(x=cancer_type, y=expression, group=interaction(cancer_type, gene))) +
    geom_boxplot(aes(fill=gene), outlier.size=0.1, position=position_dodge(0.7), width=0.5, lwd=0.1) +
    fill_palette(PAL_GENE) +
    theme_pubr() +
    geom_text(aes(label=p.signif, y=9, group=gene), tests, size=2, 
              family=FONT_FAMILY, position=position_dodge(0.7)) +
    labs(x="Cancer Type", y="log2(Norm. Count + 1)", fill="Gene")

x = x %>% mutate(cancer_type="PANCAN")
tests = compare_means(
    formula=expression~gene, data=x, 
    method=TEST_METHOD, group.by="cancer_type") %>%
    mutate(gene=group2)
plts[["expression-boxplot-pancan"]] = x %>%
    ggplot(aes(x=cancer_type, y=expression, group=interaction(cancer_type, gene))) +
    geom_boxplot(aes(fill=gene), outlier.size=0.1, position=position_dodge(0.7), width=0.5, lwd=0.1) +
    fill_palette(PAL_GENE) +
    theme_pubr() +
    geom_text(aes(label=p.signif, y=9, group=gene), tests, size=2, 
              family=FONT_FAMILY, position=position_dodge(0.7)) +
    labs(x="Cancer Type", y="log2(Norm. Count + 1)", fill="Gene")


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
    save_plot(filename, plt, base_width=width, base_height=height, dpi=dpi, units='cm')
}

## plots
dir.create(figs_dir, recursive=TRUE)
save_plt(plts, "coexpression-scatter", ".pdf", figs_dir, width=10, height=12)
save_plt(plts, "expression-boxplot-bycancer", ".pdf", figs_dir, width=12, height=7)
save_plt(plts, "expression-boxplot-pancan", ".pdf", figs_dir, width=2.25, height=7)

# figdata
filename = file.path(RESULTS_DIR,"files","figdata-TTLL11_vs_TTLL13-genexpr.xlsx")
X %>%
    pivot_longer(c(TTLL11,TTLL13), names_to="gene", values_to="expression") %>%
    write_xlsx(filename)


print("Done!")