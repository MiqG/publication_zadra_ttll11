
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
corrs = lapply(c("Primary Tumor","Solid Tissue Normal"), function(sample_type_oi){
    samples_oi = metadata %>% filter(sample_type==sample_type_oi) %>% pull(sample) %>% unique()
    common_samples = intersect(colnames(genexpr), samples_oi)

    corr = cor(genexpr[,common_samples] %>% t(), 
               genexpr_oi[,common_samples] %>% t(), 
               method="spearman") %>%
        as.data.frame() %>%
        mutate(symbol = genexpr[["sample"]]) %>%
        drop_na() %>%
        pivot_longer(-symbol, names_to="gene", values_to="correlation") %>% 
        filter(symbol != gene) %>% # do not consider self-correlations
        mutate(sample_type=sample_type_oi)
    
    return(corr)
})
corrs = do.call(rbind,corrs)

## Primary Tumors
X = corrs %>% filter(sample_type=="Primary Tumor")

plts[["coexpression_tcga-TTLLs_vs_all-violin-pt"]] = X %>% 
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

X = corrs %>% filter(sample_type=="Solid Tissue Normal")

plts[["coexpression_tcga-TTLLs_vs_all-violin-stn"]] = X %>% 
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

X = corrs
plts[["coexpression_tcga-TTLLs_vs_all-scatter-pt_vs_stn"]] = X %>% 
    pivot_wider(id_cols=c("symbol","gene"), names_from="sample_type", 
                values_from="correlation", values_fn = mean) %>%
    drop_na() %>%
    ggplot(aes(x=`Primary Tumor`, y=`Solid Tissue Normal`)) +
    geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
    facet_wrap(~gene) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    theme_pubr()

# gene set enrichment analysis for each TTLL
enrichments = list()

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

## Primary Tumor
X = corrs %>% filter(sample_type=="Primary Tumor")
query = sapply(GENES_OI, function(gene_oi){
    X %>% 
        filter(gene %in% gene_oi) %>% 
        distinct(symbol,correlation) %>% 
        arrange(-correlation) %>% 
        deframe()
}, simplify=FALSE)

### GO
enrichments[["PT"]][["GO"]] = sapply(GENES_OI, function(gene_oi){
    gseGO(query[[gene_oi]], ont="BP", OrgDb=ORGDB, keyType="SYMBOL")
}, simplify=FALSE)
enrichments[["PT"]][["GO"]] = get_enrichment_result(enrichments[["PT"]][["GO"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)

### CHEA TF targets
enrichments[["PT"]][["tf_targets"]] = sapply(GENES_OI, function(gene_oi){
    GSEA(
        query[[gene_oi]], maxGSSize=1e4, 
        TERM2GENE=ontologies[["tf_targets"]] %>% 
            group_by(term) %>% 
            filter(any(gene %in% GENES_OI)) %>%
            ungroup()
    )
}, simplify=FALSE)
enrichments[["PT"]][["tf_targets"]] = get_enrichment_result(enrichments[["PT"]][["tf_targets"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)


## Solid Tissue Normal
X = corrs %>% filter(sample_type=="Solid Tissue Normal")
query = sapply(GENES_OI, function(gene_oi){
    X %>% 
        filter(gene %in% gene_oi) %>% 
        distinct(symbol,correlation) %>% 
        arrange(-correlation) %>% 
        deframe()
}, simplify=FALSE)

### GO
enrichments[["STN"]][["GO"]] = sapply(GENES_OI, function(gene_oi){
    gseGO(query[[gene_oi]], ont="BP", OrgDb=ORGDB, keyType="SYMBOL")
}, simplify=FALSE)
enrichments[["STN"]][["GO"]] = get_enrichment_result(enrichments[["STN"]][["GO"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)

### CHEA TF targets
enrichments[["STN"]][["tf_targets"]] = sapply(GENES_OI, function(gene_oi){
    GSEA(
        query[[gene_oi]], maxGSSize=1e4, 
        TERM2GENE=ontologies[["tf_targets"]] %>% 
            group_by(term) %>% 
            filter(any(gene %in% GENES_OI)) %>%
            ungroup()
    )
}, simplify=FALSE)
enrichments[["STN"]][["tf_targets"]] = get_enrichment_result(enrichments[["STN"]][["tf_targets"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)


# visualize enrichments
## Primary Tumor
### GO
terms_oi = enrichments[["PT"]][["GO"]] %>% 
    group_by(Cluster) %>%
    slice_min(NES, n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES-pt"]] = enrichments[["PT"]][["GO"]] %>% 
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

terms_oi = enrichments[["PT"]][["GO"]] %>% 
    group_by(Cluster) %>%
    slice_max(NES, n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES-pt"]] = enrichments[["PT"]][["GO"]] %>% 
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

#### CHEA TFs
terms_oi = enrichments[["PT"]][["tf_targets"]] %>% 
    group_by(Cluster) %>%
    slice_max(NES, n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES-pt"]] = enrichments[["PT"]][["tf_targets"]] %>% 
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)


## Solid Tissue Normal
### GO
terms_oi = enrichments[["STN"]][["GO"]] %>% 
    group_by(Cluster) %>%
    slice_min(NES, n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES-stn"]] = enrichments[["STN"]][["GO"]] %>% 
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

terms_oi = enrichments[["STN"]][["GO"]] %>% 
    group_by(Cluster) %>%
    slice_max(NES, n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES-stn"]] = enrichments[["STN"]][["GO"]] %>% 
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

#### CHEA TFs
terms_oi = enrichments[["STN"]][["tf_targets"]] %>% 
    group_by(Cluster) %>%
    slice_max(NES, n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES-stn"]] = enrichments[["STN"]][["tf_targets"]] %>% 
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="GeneRatio", color="p.adjust") +
    scale_size(range=c(0.5,3)) + 
    scale_color_continuous(
        low=PAL_FDR_LIGHT, high=PAL_FDR_DARK, 
        name="FDR", guide=guide_colorbar(reverse=TRUE)) +
    theme_pubr(x.text.angle = 70)

## PT vs STN
### GO
X = enrichments[["STN"]][["GO"]] %>%
    distinct(Description, NES, Cluster) %>%
    left_join(
        enrichments[["PT"]][["GO"]] %>%
        distinct(Description, NES, Cluster),
        by = c("Description","Cluster"),
        suffix = c("_stn","_pt")
    ) %>%
    drop_na()

plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment-pt_vs_stn"]] = X %>%
    ggplot(aes(x=NES_pt, y=NES_stn)) +
    geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
    facet_wrap(~Cluster) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    theme_pubr() +
    labs(x="NES Primary Tumor", y="NES Solid Tissue Normal")


### CHEA TFs
X = enrichments[["STN"]][["tf_targets"]] %>%
    distinct(Description, NES, Cluster) %>%
    left_join(
        enrichments[["PT"]][["tf_targets"]] %>%
        distinct(Description, NES, Cluster),
        by = c("Description","Cluster"),
        suffix = c("_stn","_pt")
    ) %>%
    drop_na()

plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-pt_vs_stn"]] = X %>%
    ggplot(aes(x=NES_pt, y=NES_stn)) +
    geom_scattermore(pixels=c(1000,1000), pointsize=8, alpha=0.5) +
    facet_wrap(~Cluster) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    theme_pubr() +
    labs(x="NES Primary Tumor", y="NES Solid Tissue Normal")


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
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-violin-pt", ".pdf", figs_dir, width=15, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-violin-stn", ".pdf", figs_dir, width=15, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-scatter-pt_vs_stn", ".pdf", figs_dir, width=15, height=10)

save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES-pt", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES-pt", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES-pt", ".pdf", figs_dir, width=10, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES-stn", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES-stn", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES-stn", ".pdf", figs_dir, width=10, height=10)

save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment-pt_vs_stn", ".pdf", figs_dir, width=10, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-pt_vs_stn", ".pdf", figs_dir, width=10, height=10)

print("Done!")