
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
RAW_DIR = file.path(ROOT,'data','raw')
PREP_DIR = file.path(ROOT,'data','prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
ORGDB = org.Hs.eg.db
ORDER_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP',
             'THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
TEST_METHOD = "wilcox.test"
GENE_OI = "TTLL11"
GENES_OI = c('TTLL1', 'TTLL2', 'TTLL4', 'TTLL5', 'TTLL6', 'TTLL7', 'TTLL9', 'TTLL11', 'TTLL13')

# formatting
FONT_SIZE = 7 # pt
FONT_FAMILY = 'Helvetica'
PAL_SAMPLE_TYPE = rev(get_palette("npg",2))
PAL_FDR_DARK = "#005AB5"
PAL_FDR_LIGHT = "#DC3220"

# inputs
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
genexpr_file = file.path(PREP_DIR,'genexpr_TCGA.tsv.gz')
ontology_chea_file = file.path(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz")

# outputs
figs_dir = file.path(RESULTS_DIR,"figures","coexpression_tcga")

# load
metadata = read_tsv(phenotype_file)
genexpr = read_tsv(genexpr_file)
ontologies = list(
    "tf_targets" = read.gmt(ontology_chea_file)
)

# prep
metadata = metadata %>%
    filter(cancer_type %in% ORDER_OI)

genexpr_oi = genexpr %>% 
    dplyr::rename(gene=sample) %>% 
    filter(gene %in% GENES_OI) %>% 
    column_to_rownames("gene")

plts = list()

# does the expresson of TTLL11 correlate with expression of any other gene?

tumor_samples = metadata %>% filter(sample_type=="Primary Tumor") %>% pull(sample) %>% unique()
common_samples = intersect(colnames(genexpr), tumor_samples)

corr = cor(genexpr[,common_samples] %>% t(), 
           genexpr_oi[,common_samples] %>% t(), 
           method="spearman") %>%
    as.data.frame() %>%
    mutate(symbol = genexpr[["sample"]]) %>%
    drop_na() %>%
    pivot_longer(-symbol, names_to="gene", values_to="correlation") %>% 
    filter(symbol != gene) # do not consider self-correlations

X = corr 

plts[["coexpression_tcga-TTLLs_vs_all-violin"]] = X %>% 
    ggviolin(x="gene", y="correlation", color=NA, fill="orange") + 
    geom_boxplot(width=0.1, outlier.size=0.1) + 
    labs(x="Gene", y="Correlation", subtitle=sprintf("n=%s",nrow(X))) + 
    geom_text_repel(aes(label=symbol), 
                    X %>% group_by(gene) %>% slice_max(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY) + 
    geom_text_repel(aes(label=symbol), 
                    X %>% group_by(gene) %>% slice_min(correlation, n=10), 
                    max.overlaps=50, segment.size=0.1,
                    size=2, family=FONT_FAMILY)

# gene set enrichment analysis for each TTLL
query = sapply(GENES_OI, function(gene_oi){
    X %>% 
        filter(gene %in% gene_oi) %>% 
        distinct(symbol,correlation) %>% 
        arrange(-correlation) %>% 
        deframe()
}, simplify=FALSE)

get_enrichment_result = function(enrich_list, thresh=0.05){
    ## groups are extracted from names
    groups = names(enrich_list)
    result = lapply(groups, function(group){
        res = enrich_list[[group]]
        if(nrow(res@result)>0){
            res = res@result
            res$Cluster = group
        }else{
            res = NULL
        }
        return(res)
    })
    result[sapply(result, is.null)] = NULL
    result = do.call(rbind,result)
    
    ## filter by p.adjusted
    result = result %>% filter(p.adjust<thresh)
    
    return(result)
}

enrichments = list()
## GO
enrichments[["GO"]] = sapply(GENES_OI, function(gene_oi){
    gseGO(query[[gene_oi]], ont="BP", OrgDb=ORGDB, keyType="SYMBOL")
}, simplify=FALSE)
enrichments[["GO"]] = get_enrichment_result(enrichments[["GO"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)

plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES"]] = enrichments[["GO"]] %>% 
    group_by(Cluster) %>%
    slice_min(NES, n=5) %>%
    ungroup() %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES"]] = enrichments[["GO"]] %>% 
    group_by(Cluster) %>%
    slice_max(NES, n=5) %>%
    ungroup() %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

## CHEA TF targets
enrichments[["tf_targets"]] = sapply(GENES_OI, function(gene_oi){
    GSEA(
        query[[gene_oi]], maxGSSize=1e4, 
        TERM2GENE=ontologies[["tf_targets"]] %>% 
            group_by(term) %>% 
            filter(any(gene %in% GENES_OI)) %>%
            ungroup()
    )
}, simplify=FALSE)
enrichments[["tf_targets"]] = get_enrichment_result(enrichments[["tf_targets"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)

plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES"]] = enrichments[["tf_targets"]] %>% 
    group_by(Cluster) %>%
    slice_max(NES, n=5) %>% # positively coexpressed should be co-regulated
    ungroup() %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)


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
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-violin", ".pdf", figs_dir, width=15, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES", ".pdf", figs_dir, width=10, height=10)
