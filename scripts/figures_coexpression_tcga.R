
require(tidyverse)
require(ggpubr)
require(cowplot)
require(writexl)
require(scattermore)
require(ggrepel)
require(clusterProfiler)
require(org.Hs.eg.db)
require(enrichplot)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
ORGDB = org.Hs.eg.db
ORDER_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP',
             'THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
TEST_METHOD = "wilcox.test"
GENE_OI = "TTLL11"

# formatting
FONT_SIZE = 7 # pt
FONT_FAMILY = 'Helvetica'
PAL_SAMPLE_TYPE = rev(get_palette("npg",2))
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(PREP_DIR,'genexpr_TCGA.tsv.gz')

# outputs
figs_dir = file.path(RESULTS_DIR,"figures","coexpression_tcga")

# load
metadata = read_tsv(phenotype_file)
genexpr = read_tsv(genexpr_file)

# prep
metadata = metadata %>%
    filter(cancer_type %in% ORDER_OI)

genexpr_oi = genexpr %>% 
    rename(gene=sample) %>% 
    filter(gene %in% GENE_OI) %>% 
    pivot_longer(-gene, names_to="sample", values_to="expression") %>%
    dplyr::select("sample", "expression") %>% 
    deframe()

plts = list()

# does the expresson of TTLL11 correlate with expression of any other gene?

tumor_samples = metadata %>% filter(sample_type=="Primary Tumor") %>% pull(sample) %>% unique()
common_samples = intersect(colnames(genexpr), tumor_samples)

corr = cor(genexpr[,common_samples] %>% t(), 
           genexpr_oi[common_samples], 
           method="spearman", use="pairwise.complete.obs") %>%
    as.vector()

corr = data.frame(
    symbol = genexpr[["sample"]],
    correlation = corr
)

X = corr %>% 
    mutate(gene = "TTLL11") %>%
    filter(symbol != "TTLL11")

plts[["coexpression_tcga-TTLL11_vs_all-violin"]] = X %>% 
    ggviolin(x="gene", y="correlation", color=NA, fill="orange") + 
    geom_boxplot(width=0.1, outlier.size=0.1) + 
    labs(x="Gene", y="Correlation", subtitle=sprintf("n=%s",nrow(X))) + 
    geom_text_repel(aes(label=symbol), 
                    X %>% slice_max(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY) + 
    geom_text_repel(aes(label=symbol), 
                    X %>% slice_min(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY)

# gene set enrichment analysis
query = X %>% drop_na() %>% distinct(symbol,correlation) %>% arrange(-correlation) %>% deframe()
result = gseGO(query, ont="BP", OrgDb=ORGDB, keyType="SYMBOL")
plts[["coexpression_tcga-TTLL11_vs_all-enrichment_dotplot"]] = result %>% 
    dotplot() + 
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr()

# enrichment with transcription factor term
genes_oi = result %>% 
    as.data.frame() %>% 
    filter(str_detect(Description, "transcription factor")) %>% 
    slice_min(setSize, n=1) %>%
    pull(core_enrichment) %>% 
    str_split("/") %>% 
    unlist()

x = X %>% 
    filter(symbol %in% genes_oi) %>% 
    arrange(-correlation) 
plts[["coexpression_tcga-TTLL11_vs_all-enrichment_TFs"]] = x %>% 
    ggviolin(x="gene", y="correlation", color=NA, fill="orange") + 
    geom_boxplot(width=0.1, outlier.size=0.1) + 
    labs(x="Gene", y="Correlation", subtitle=sprintf("n=%s", length(genes_oi))) + 
    geom_text_repel(aes(label=symbol), 
                    x %>% slice_max(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY) + 
    geom_text_repel(aes(label=symbol), 
                    x %>% slice_min(correlation, n=10), 
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
save_plt(plts, "coexpression_tcga-TTLL11_vs_all-violin", ".pdf", figs_dir, width=7, height=7)
save_plt(plts, "coexpression_tcga-TTLL11_vs_all-enrichment_dotplot", ".pdf", figs_dir, width=12, height=7)
save_plt(plts, "coexpression_tcga-TTLL11_vs_all-enrichment_TFs", ".pdf", figs_dir, width=7, height=7)
