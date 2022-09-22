#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compare gene expression of TTLLs in TCGA.
#

require(tidyverse)
require(cowplot)
require(ggpubr)
require(writexl)
require(latex2exp)
require(scattermore)
require(extrafont)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
TCGA_DIR = file.path(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = file.path(DATA_DIR,'prep')
RESULTS_DIR = file.path(ROOT,'results')

# variables
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP',
               'THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')

THRESH_LOW = 0.10
THRESH_HIGH = 0.90

# THRESH_LOW_CYCLIN = 5
# THRESH_HIGH_CYCLIN = 10
# THRESH_LOW_CDC25 = 3.5
# THRESH_HIGH_CDC25 = 8.5

TEST_METHOD = 'wilcox.test'

# formatting
PAL_SAMPLE_TYPE = setNames(get_palette("npg",2), c("Primary Tumor", "Solid Tissue Normal"))
FONT_SIZE = 7 # pt
FONT_FAMILY = 'Arial'
PAL_STATUS = setNames(
    c("orange","darkred","black",'#4DBBD5FF'),
    #c("#E9B872","#A63D40","#151515"),
    c("CCNE1_low & CDC25A_low", "CCNE1_high or CDC25A_high", "N.R.","STN")
)

# inputs
genexpr_ttlls_file = file.path(PREP_DIR,'genexpr_TTLLs.tsv')
genexpr_file = file.path(PREP_DIR,'genexpr_TCGA.tsv.gz')
phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')

# outputs
figs_dir = file.path(RESULTS_DIR,"figures","differential_genexpr_by_oncogene_status")

# load data
genexpr_ttlls = read_tsv(genexpr_ttlls_file)
genexpr = read_tsv(genexpr_file)
metadata = read_tsv(phenotype_file)

# select genes of interest and prepare
X = genexpr %>% 
    filter(sample %in% c("CCNE1","CDC25A")) %>%
    column_to_rownames("sample") %>% 
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(metadata, by="sample") %>% 
    mutate(sample_type = factor(sample_type, levels = SAMPLE_TYPES_OI),
           cancer_type = factor(cancer_type, levels = CANCERS_OI)) %>%
    drop_na()

plts = list()

##### How does the expression of CyclinE and Cdc25 change across cancer types? #####
# classify samples based on high or low expression of the genes
X = X %>%
    group_by(cancer_type) %>%
    mutate(
        status_cyclin = case_when(
            CCNE1 < quantile(CCNE1, THRESH_LOW) ~ "CCNE1_low",
            CCNE1 > quantile(CCNE1, THRESH_HIGH) ~ "CCNE1_high",
            CCNE1 >= quantile(CCNE1, THRESH_LOW) & CCNE1 <= quantile(CCNE1, THRESH_HIGH) ~ "CCNE1_regular"
        ),
        status_cdc25 = case_when(
            CDC25A < quantile(CDC25A, THRESH_LOW) ~ "CDC25A_low",
            CDC25A > quantile(CDC25A, THRESH_HIGH) ~ "CDC25A_high",
            CDC25A >= quantile(CDC25A, THRESH_LOW) & CDC25A <= quantile(CDC25A, THRESH_HIGH) ~ "CDC25A_regular"
        ),
        status_combined = paste0(status_cyclin, " & ", status_cdc25),
    ) %>%
    ungroup() %>%
    mutate(
        status_combined_clean = case_when(
            status_cdc25 == "CDC25A_high" | status_cdc25 == "CCNE1_high" ~ "CCNE1_high or CDC25A_high",
            status_combined == "CCNE1_low & CDC25A_low" ~ "CCNE1_low & CDC25A_low",
            TRUE ~ "N.R."
        ),
        status_combined_clean = ifelse(sample_type == "Solid Tissue Normal", "STN", status_combined_clean)
    )

# viz
status_oi = c("CCNE1_low & CDC25A_low", "CCNE1_high or CDC25A_high", "N.R.")
plts[["coexpression_cyclinE_cdc25"]] = X %>% 
    filter(sample_type=="Primary Tumor") %>%
    ggplot(aes(x=CCNE1, y=CDC25A)) + 
    geom_scattermore(aes(color=status_combined_clean), pixels=c(1000,1000), pointsize=10, alpha=0.5) +
    stat_cor(method="spearman", size=1.5, family=FONT_FAMILY) +
    theme_pubr() +
    facet_wrap(~cancer_type, scales="free") +
    theme(aspect.ratio=1, strip.text.x = element_text(size=6, family=FONT_FAMILY)) +
    color_palette(PAL_STATUS[status_oi]) +
    labs(x="log2(Norm. Count CCNE1 + 1)", y="log2(Norm. Count CDC25A + 1)", color="Gene Status")

##### Run Differential analysis for all TTLLs across cancers considering oncogenes status #####
# which cancers have enough samples?
comparisons_oi = X %>% 
    filter(sample_type=="Primary Tumor") %>% 
    count(cancer_type, status_combined_clean) %>% 
    filter(!(status_combined_clean=="N.R.")) #& n>20)

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

# run analysis with TTLLs (TO SAVE)
df = comparisons_oi %>%
    left_join(X, by=c("cancer_type","status_combined_clean")) %>%
    bind_rows(X %>% filter(cancer_type%in%comparisons_oi[["cancer_type"]] & status_combined_clean=="STN")) %>%
    distinct(cancer_type, status_combined_clean, sample, sample_type) %>%
    left_join(genexpr_ttlls, by=c("sample","cancer_type","sample_type"))

cancers_oi = df %>% pull(cancer_type) %>% unique()
result = lapply(cancers_oi, function(cancer_type_oi){
    print(cancer_type_oi)
    
    status_oi = df %>% 
        filter(cancer_type %in% cancer_type_oi) %>% 
        pull(status_combined_clean) %>% 
        unique() %>%
        setdiff("STN")
    
    result = lapply(status_oi, function(status_combined_oi){
        tmp = df %>% 
            filter(cancer_type%in%cancer_type_oi & 
                   (status_combined_clean %in% c(status_combined_oi,"STN")))
        
        result = run_wilcox(tmp) %>% 
            mutate(cancer_type = cancer_type_oi,
                   status_combined_clean = status_combined_oi)
        
        return(result)
    })
    result = do.call(rbind, result)

}) 
result = do.call(rbind, result)

##### Plot Gene Expression Distributions of TTLL11 across cancers #####
status_order = c(
    'STN',
    'CCNE1_high & CDC25A_high',
    'CCNE1_high & CDC25A_low',
    'CCNE1_low & CDC25A_high',
    'CCNE1_high & CDC25A_regular',
    'CCNE1_regular & CDC25A_high',
    'CCNE1_regular & CDC25A_low',
    'CCNE1_low & CDC25A_regular',
    'CCNE1_low & CDC25A_low'
)

status_order = c(
    'STN',
    'CCNE1_high or CDC25A_high',
    'CCNE1_low & CDC25A_low'
)

plts[["differential_expression_by_status-TTLL11"]] = df %>%
    filter(gene == "TTLL11" & status_combined_clean%in%status_order) %>%
    mutate(status_combined_clean = factor(status_combined_clean, levels=status_order),
           cancer_type = factor(cancer_type, levels=CANCERS_OI)) %>%
    ggplot(aes(x=status_combined_clean, y=expression)) +
    geom_boxplot(aes(fill=sample_type), outlier.size=0.1, lwd=0.1) +
    fill_palette(PAL_SAMPLE_TYPE) +
    theme_pubr(x.text.angle = 70) +
    facet_wrap(~cancer_type) +
    theme(strip.text.x = element_text(size=6, family=FONT_FAMILY),
          aspect.ratio = 1) +
    stat_compare_means(aes(label=..p.signif..),
                       method=TEST_METHOD, ref.group="STN", size=2, family=FONT_FAMILY) + 
    labs(x="Oncogene Status", y=TeX('$log_2(Norm. Count + 1)$'), fill="Sample Type")


comparisons = list(
    c("STN","CCNE1_high or CDC25A_high"),
    c("STN","CCNE1_low & CDC25A_low")
)
x = df %>%
    filter(gene == "TTLL11" & status_combined_clean%in%status_order) %>%
    mutate(status_combined_clean = factor(status_combined_clean, levels=status_order),
           cancer_type = factor(cancer_type, levels=CANCERS_OI))
tests = compare_means(
    formula=expression~status_combined_clean, data=x, 
    method=TEST_METHOD, ref.group="STN", 
    group.by="cancer_type", comparisons=comparisons) %>%
    mutate(status_combined_clean=group2)

plts[["differential_expression_by_status_and_cancer-TTLL11"]] = x %>%
    ggplot(aes(x=cancer_type, y=expression, group=interaction(cancer_type,status_combined_clean))) +
    geom_boxplot(aes(fill=status_combined_clean), outlier.size=0.1, position=position_dodge(0.7), width=0.5, lwd=0.1) +
    fill_palette(PAL_STATUS[status_order]) +
    theme_pubr() +
    geom_text(aes(label=p.signif, y=9, group=status_combined_clean), tests, size=2, 
              family=FONT_FAMILY, position=position_dodge(0.7)) +
    labs(x="Cancer Type", y="log2(Norm. Count + 1)", fill="Sample Type")

x = x %>% mutate(cancer_type="PANCAN")
tests = compare_means(
    formula=expression~status_combined_clean, data=x, 
    method=TEST_METHOD, ref.group="STN", 
    group.by="cancer_type", comparisons=comparisons) %>%
    mutate(status_combined_clean=group2)

plts[["differential_expression_by_status_and_pancan-TTLL11"]] = x %>%
    ggplot(aes(x=cancer_type, y=expression, group=interaction(cancer_type,status_combined_clean))) +
    geom_boxplot(aes(fill=status_combined_clean), outlier.size=0.1, position=position_dodge(0.7), width=0.5, lwd=0.1) +
    fill_palette(PAL_STATUS[status_order]) +
    theme_pubr() +
    geom_text(aes(label=p.signif, y=9, group=status_combined_clean), tests, size=2, 
              family=FONT_FAMILY, position=position_dodge(0.7)) +
    labs(x="Cancer Type", y="log2(Norm. Count + 1)", fill="Sample Type")

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
save_plt(plts, "coexpression_cyclinE_cdc25", ".pdf", figs_dir, width=10, height=12)
save_plt(plts, "differential_expression_by_status-TTLL11", ".pdf", figs_dir, width=15, height=18)
save_plt(plts, "differential_expression_by_status_and_cancer-TTLL11", ".pdf", figs_dir, width=12, height=7)
save_plt(plts, "differential_expression_by_status_and_pancan-TTLL11", ".pdf", figs_dir, width=2.25, height=7)

## figdata
filename = file.path(RESULTS_DIR,"files","figdata-ccne1_cdc25a_tcga-genexpr.xlsx")
write_xlsx(X, filename)
filename = file.path(RESULTS_DIR,"files","figdata-ccne1_cdc25a_tcga-diffgenexpr_TTLL11.xlsx")
result %>% filter(gene=="TTLL11") %>% write_xlsx(filename)


print("Done!")