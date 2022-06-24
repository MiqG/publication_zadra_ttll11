
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
REGULATORS = c("CDC25A","CCNE1")

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
    
    gc()
    
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
    res = gseGO(query[[gene_oi]], ont="BP", OrgDb=ORGDB, keyType="SYMBOL")
    gc()
    return(res)
}, simplify=FALSE)
enrichments[["PT"]][["GO"]] = get_enrichment_result(enrichments[["PT"]][["GO"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)
gc()

### CHEA TF targets
enrichments[["PT"]][["tf_targets"]] = sapply(GENES_OI, function(gene_oi){
    res = GSEA(
        query[[gene_oi]], maxGSSize=1e4, 
        TERM2GENE=ontologies[["tf_targets"]] %>% 
            group_by(term) %>% 
            filter(any(gene %in% GENES_OI)) %>%
            ungroup()
    )
    gc()
    return(res)
}, simplify=FALSE)
enrichments[["PT"]][["tf_targets"]] = get_enrichment_result(enrichments[["PT"]][["tf_targets"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)
gc()

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
    res = gseGO(query[[gene_oi]], ont="BP", OrgDb=ORGDB, keyType="SYMBOL")
    gc()
    return(res)
}, simplify=FALSE)
enrichments[["STN"]][["GO"]] = get_enrichment_result(enrichments[["STN"]][["GO"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)
gc()

### CHEA TF targets
enrichments[["STN"]][["tf_targets"]] = sapply(GENES_OI, function(gene_oi){
    res = GSEA(
        query[[gene_oi]], maxGSSize=1e4, 
        TERM2GENE=ontologies[["tf_targets"]] %>% 
            group_by(term) %>% 
            filter(any(gene %in% GENES_OI)) %>%
            ungroup()
    )
    gc()
    return(res)
}, simplify=FALSE)
enrichments[["STN"]][["tf_targets"]] = get_enrichment_result(enrichments[["STN"]][["tf_targets"]]) %>%
    mutate(Count = str_count(core_enrichment, "/")+1, 
           GeneRatio=Count/setSize)
gc()

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
    distinct(Description, NES, Cluster, GeneRatio, core_enrichment) %>%
    left_join(
        enrichments[["PT"]][["GO"]] %>%
        distinct(Description, NES, Cluster, GeneRatio, core_enrichment),
        by = c("Description","Cluster"),
        suffix = c("_stn","_pt")
    ) %>%
    mutate(NES_diff = NES_pt - NES_stn,
           abs_NES_diff = abs(NES_diff)) %>%
    drop_na()

plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment-pt_vs_stn"]] = X %>%
    ggplot(aes(x=NES_pt, y=NES_stn)) +
    geom_scattermore(pixels=c(1000,1000), pointsize=4, alpha=0.5) +
    facet_wrap(~Cluster) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    theme_pubr() +
    labs(x="NES Primary Tumor", y="NES Solid Tissue Normal")

terms_oi = X %>% 
    group_by(Cluster) %>%
    slice_max(abs(NES_diff), n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment-pt_vs_stn_diff"]] = X %>%
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="abs_NES_diff", color="NES_diff") +
    scale_size("|Delta NES|", range=c(0.5,3)) + 
    scale_color_continuous(low="darkgreen", high="orange", name="Delta NES") +
    theme_pubr(x.text.angle = 70)

#### which gene sets contain drivers of lowering TTLL11
x = X %>% 
    filter(NES_diff<(-1) & Cluster=="TTLL11") %>% # 103 gene sets
    filter(
        str_detect(string=core_enrichment_stn, pattern=paste(REGULATORS, collapse="|")) | 
        str_detect(string=core_enrichment_pt, pattern=paste(REGULATORS, collapse="|"))
    ) %>%
    rowwise() %>%
    mutate(
        matches_stn = paste(intersect(REGULATORS, unlist(strsplit(core_enrichment_stn, "/"))), 
                            collapse=";"),
        matches_pt = paste(intersect(REGULATORS, unlist(strsplit(core_enrichment_pt, "/"))), 
                           collapse=";")
    ) %>%
    ungroup() # 19 gene sets
# 19 out of 103 differential gene sets contain regulators
# only PT samples contain regulators as core enrichment

plts[["coexpression_tcga-TTLLs_vs_all-go_enrichment-regulators_enriched"]] = x %>%
    separate_rows(matches_stn, matches_pt) %>%
    distinct(Description, matches_stn, matches_pt) %>%
    pivot_longer(-Description, names_to="sample_type", values_to="regulators") %>%
    mutate(regulators = ifelse(regulators == "", NA, regulators)) %>%
    drop_na(regulators) %>%
    count(sample_type, regulators, Description) %>%
    ggbarplot(x="regulators", y="n", fill="Description", palette=get_palette("Paired",19), color=NA) +
    facet_wrap(~sample_type) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    labs(x="Regulator TTLL11", y="Count in Gene Sets")


### CHEA TFs
X = enrichments[["STN"]][["tf_targets"]] %>%
    distinct(Description, NES, Cluster, GeneRatio, core_enrichment) %>%
    left_join(
        enrichments[["PT"]][["tf_targets"]] %>%
        distinct(Description, NES, Cluster, GeneRatio, core_enrichment),
        by = c("Description","Cluster"),
        suffix = c("_stn","_pt")
    ) %>%
    mutate(NES_diff = NES_pt - NES_stn,
           abs_NES_diff = abs(NES_diff)) %>%
    drop_na()

plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-pt_vs_stn"]] = X %>%
    ggplot(aes(x=NES_pt, y=NES_stn)) +
    geom_scattermore(pixels=c(1000,1000), pointsize=8, alpha=0.5) +
    facet_wrap(~Cluster) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    stat_cor(method="spearman", size=2, family=FONT_FAMILY) +
    theme_pubr() +
    labs(x="NES Primary Tumor", y="NES Solid Tissue Normal")

terms_oi = X %>% 
    group_by(Cluster) %>%
    slice_max(abs(NES_diff), n=5) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-pt_vs_stn_diff"]] = X %>%
    filter(Description %in% terms_oi) %>%
    arrange(desc(Cluster), Description) %>%
    mutate(Description = factor(Description, levels=unique(Description))) %>%
    ggscatter(x="Cluster", y="Description", size="abs_NES_diff", color="NES_diff") +
    scale_size("|Delta NES|", range=c(0.5,3)) + 
    scale_color_continuous(low="darkgreen", high="orange", name="Delta NES") +
    theme_pubr(x.text.angle = 70)

#### which gene sets contain drivers of lowering TTLL11
x = X %>% 
    filter(NES_diff<(-1) & Cluster=="TTLL11") %>% # 3 gene sets
    filter(
        str_detect(string=core_enrichment_stn, pattern=paste(REGULATORS, collapse="|")) | 
        str_detect(string=core_enrichment_pt, pattern=paste(REGULATORS, collapse="|"))
    ) %>%
    rowwise() %>%
    mutate(
        matches_stn = paste(intersect(REGULATORS, unlist(strsplit(core_enrichment_stn, "/"))), 
                            collapse=";"),
        matches_pt = paste(intersect(REGULATORS, unlist(strsplit(core_enrichment_pt, "/"))), 
                           collapse=";")
    ) %>%
    ungroup() # 3 gene sets
# 3 out of 3 differential gene sets contain regulators
# only PT samples contain regulators as core enrichment

plts[["coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-regulators_enriched"]] = x %>%
    separate_rows(matches_stn, matches_pt) %>%
    distinct(Description, matches_stn, matches_pt) %>%
    pivot_longer(-Description, names_to="sample_type", values_to="regulators") %>%
    mutate(regulators = ifelse(regulators == "", NA, regulators)) %>%
    drop_na(regulators) %>%
    count(sample_type, regulators, Description) %>%
    ggbarplot(x="regulators", y="n", fill="Description", palette=get_palette("Dark2",3), color=NA) +
    facet_wrap(~sample_type) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY)) + 
    labs(x="Regulator TTLL11", y="Count in Gene Sets")

# coexpression of regulators with TTLLs
X = corrs
plts[["coexpression_tcga-TTLLs_vs_all-regulators"]] = X %>% 
    filter(symbol%in%REGULATORS & gene=="TTLL11") %>% 
    ggstripchart(x="gene", y="correlation", color="sample_type", 
                 shape="symbol", palette=PAL_SAMPLE_TYPE) + 
    labs(x="TTLL Genes", y="Spearman Correlation", 
         color="Sample Type", shape="Regulator of TTLL11")

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

# separate enrichments
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES-pt", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES-pt", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES-pt", ".pdf", figs_dir, width=10, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-lowNES-stn", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment_dotplot-highNES-stn", ".pdf", figs_dir, width=15, height=12)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment_dotplot-highNES-stn", ".pdf", figs_dir, width=10, height=10)

# compare enrichments
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment-pt_vs_stn", ".pdf", figs_dir, width=10, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-go_enrichment-pt_vs_stn_diff", ".pdf", figs_dir, width=15, height=15)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-pt_vs_stn", ".pdf", figs_dir, width=10, height=10)
save_plt(plts, "coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-pt_vs_stn_diff", ".pdf", figs_dir, width=8, height=10)

# TTLL11 regulators in differentially enriched coexpression gene sets
save_plt(plts, 'coexpression_tcga-TTLLs_vs_all-go_enrichment-regulators_enriched', ".pdf", figs_dir, width=4, height=8)
save_plt(plts, 'coexpression_tcga-TTLLs_vs_all-tf_targets_enrichment-regulators_enriched', ".pdf", figs_dir, width=4, height=6)

save_plt(plts, "coexpression_tcga-TTLLs_vs_all-regulators", ".pdf", figs_dir, width=5, height=5)


print("Done!")