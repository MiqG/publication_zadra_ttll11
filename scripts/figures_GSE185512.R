# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compare gene expression of TTLLs in TCGA.
#

require(tidyverse)
require(scattermore)
require(ggpubr)
require(ggrepel)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
RAW_DIR = file.path(DATA_DIR,"raw")
RESULTS_DIR = file.path(ROOT,'results')

genes_oi = c("ENSG00000175764")

TEST_METHOD = 'wilcox.test'


# inputs
metadata_file = file.path(RAW_DIR,"GSE185512","metadata.tsv")
counts_file = file.path(RAW_DIR,"GSE185512","genexpr_counts.tsv.gz")
dge_oe_file = file.path(RESULTS_DIR,"files","dge_overexpression-GSE185512-dox_vs_ctl-merged.tsv")
dge_t_file = file.path(RESULTS_DIR,"files","dge_overexpression-GSE185512-dox_120h_vs_48h-merged.tsv")
dge_mut_file = file.path(RESULTS_DIR,"files","dge_p53mutation-GSE185512-rpe1_120h_mut_vs_wt.tsv")

# outputs
output_figdir = file.path(RESULTS_DIR,"figures","differential_gene_expression-GSE185512")
output_figdata = file.path()

# initialize
plts = list()

# load data
metadata = read_tsv(metadata_file)
counts = read_tsv(counts_file)
gene_lengths = read_tsv("~/databases/data/GENCODE/gene_length.tsv") ###### TODO
dge_oe = read_tsv(dge_oe_file)
dge_t = read_tsv(dge_t_file)
dge_mut = read_tsv(dge_mut_file)

# prepare
metadata = metadata %>% mutate(cell_line=gsub("RPE-1","RPE1",cell_line),                               
                               cell_line=gsub("RPE1 mutTP53","RPE1-mutTP53",cell_line),
                               cell_line=gsub("RPE1-mutTP53","RPE1TP53mut",cell_line),
                               cell_line=gsub("MDA-MB231","MDAMB231",cell_line),
                               cell_line_simple=gsub("-.*","",cell_line))
dge_oe = dge_oe %>% mutate(cell_line=gsub("\\|.*","",subset_values))
dge_t = dge_t %>% mutate(cell_line=gsub("\\|.*","",subset_values))

gene_lengths = gene_lengths %>% 
    filter(!str_detect(gene,"PAR_Y")) %>% 
    mutate(gene=gsub("\\..*","",gene))

genexpr = counts %>%
    pivot_longer(-gene, names_to="sampleID", values_to="count") %>%
    left_join(gene_lengths, by="gene") %>%
    group_by(sampleID) %>%
    mutate(tpm = count / median, # median gene length
           tpm = 1e6 * tpm / sum(tpm, na.rm=TRUE),
           log2tpm = log2(tpm + 1)) %>%
    ungroup() 

genexpr_oi = genexpr %>%
    filter(gene %in% genes_oi) %>%
    left_join(metadata %>% distinct(sampleID,cell_line,cell_line_simple,treatment,timepoint,replicate), 
              by="sampleID") %>%
    arrange(cell_line)

genexpr_ctls = genexpr_oi %>%
    filter(str_detect(cell_line,"Empty")) %>%
    group_by(cell_line_simple, timepoint, treatment) %>%
    summarize(avg_ctl = mean(log2tpm, na.rm=TRUE)) %>%
    ungroup()

genexpr_oi = genexpr_oi %>%
    left_join(genexpr_ctls, by=c("cell_line_simple","timepoint","treatment")) %>%
    mutate(norm_log2tpm = log2tpm - avg_ctl)


# overview of the dataset
plts[["metadata-eda"]] = metadata %>% 
    count(cell_line, timepoint, treatment) %>% 
    ggbarplot(x="cell_line", y="n", fill="treatment", 
              palette="lancet", color=NA) + 
    facet_wrap(~timepoint) + 
    theme(strip.text.x = element_text(size=6)) + 
    labs(x="Cell Line", y="Count") + 
    coord_flip()

plts[["genexpr-eda-oe"]] = genexpr_oi %>% 
    filter(treatment!="Hydroxyurea") %>%
    mutate(treatment = factor(treatment, levels=c("NO Dox","Doxycycline"))) %>%
    ggboxplot(x="cell_line", y="log2tpm", fill="treatment", 
              palette=rev(get_palette("npg",2)), outlier.size=0.1) + 
    theme_pubr(x.text.angle = 75) +
    facet_wrap(~cell_line_simple, nrow=1, scales="free_x") +
    theme(strip.text.x = element_text(size=6)) +
    labs(x="Cell Line", y="log2(TPM + 1)", title=genes_oi)

plts[["genexpr-eda-t"]] = genexpr_oi %>% 
    filter(treatment=="Doxycycline") %>%
    mutate(timepoint = factor(timepoint, levels=c("48h","120h"))) %>%
    ggboxplot(x="cell_line", y="log2tpm", fill="timepoint", 
              palette="jco", outlier.size=0.1) + 
    facet_wrap(~treatment, ncol=1) + 
    theme_pubr(x.text.angle = 75) +
    facet_wrap(~cell_line_simple, nrow=1, scales="free_x") +
    theme(strip.text.x = element_text(size=6)) +
    labs(x="Cell Line", y="log2(TPM + 1)", title=genes_oi)

plts[["genexpr-eda-oe-norm"]] = genexpr_oi %>% 
    filter(treatment!="Hydroxyurea") %>%
    mutate(treatment = factor(treatment, levels=c("NO Dox","Doxycycline"))) %>%
    ggboxplot(x="cell_line", y="norm_log2tpm", fill="treatment", 
              palette=rev(get_palette("npg",2)), outlier.size=0.1) + 
    stat_compare_means(aes(group = treatment), 
                       method=TEST_METHOD, label = "p.signif") + 
    theme_pubr(x.text.angle = 75) +
    facet_wrap(~cell_line_simple, nrow=1, scales="free_x") +
    theme(strip.text.x = element_text(size=6)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.2) +
    labs(x="Cell Line", y="log2(TPM + 1) w.r.t. Empty", title=genes_oi)

plts[["genexpr-eda-t-norm"]] = genexpr_oi %>% 
    filter(treatment=="Doxycycline") %>%
    mutate(timepoint = factor(timepoint, levels=c("48h","120h"))) %>%
    ggboxplot(x="cell_line", y="norm_log2tpm", fill="timepoint", 
              palette="jco", outlier.size=0.1) + 
    stat_compare_means(aes(group = timepoint), 
                       method=TEST_METHOD, label = "p.signif") +
    facet_wrap(~treatment, ncol=1) + 
    theme_pubr(x.text.angle = 75) +
    facet_wrap(~cell_line_simple, nrow=1, scales="free_x") +
    theme(strip.text.x = element_text(size=6)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.2) +
    labs(x="Cell Line", y="log2(TPM + 1) w.r.t. Empty", title=genes_oi)

# effect of p53 mutation
df = dge_mut
plts[["dge-volcano-fc_vs_pvalue-mut"]] = df %>% 
    ggplot(aes(x=log2FoldChange, y=(-log10_pvalue))) + 
    geom_scattermore(pointsize=8, alpha=0.5, pixels=c(1000,1000)) + 
    geom_point(aes(x=log2FoldChange, y=(-log10_pvalue)), 
               df %>% filter(gene%in%genes_oi), 
               color="orange") +
    geom_hline(yintercept=-log10(0.05), size=0.2, linetype="dashed") + 
    labs(x="log2(FC) TP53mut vs TP53wt", y="-log10(p-value)") + 
    theme_pubr() +
    theme(strip.text.x = element_text(size=6))

plts[["dge-volcano-fc_vs_padj-mut"]] = df %>% 
    ggplot(aes(x=log2FoldChange, y=(-log10_padj))) + 
    geom_scattermore(pointsize=8, alpha=0.5, pixels=c(1000,1000)) + 
    geom_point(aes(x=log2FoldChange, y=(-log10_padj)), 
               df %>% filter(gene%in%genes_oi), 
               color="orange") +
    geom_hline(yintercept=-log10(0.1), size=0.2, linetype="dashed") + 
    labs(x="log2(FC) TP53mut vs TP53wt", y="-log10(FDR)") + 
    theme_pubr() +
    theme(strip.text.x = element_text(size=6))

# effect of overexpressing oncogenes along time
df = dge_t
plts[["dge-volcano-fc_vs_pvalue-t"]] = df %>% 
    ggplot(aes(x=log2FoldChange, y=(-log10_pvalue))) + 
    geom_scattermore(pointsize=8, alpha=0.5, pixels=c(1000,1000)) + 
    geom_point(aes(x=log2FoldChange, y=(-log10_pvalue)), 
               df %>% filter(gene%in%genes_oi), 
               color="orange") +
    geom_hline(yintercept=-log10(0.05), size=0.2, linetype="dashed") + 
    facet_wrap(~cell_line, scales="free_y") + 
    labs(x="log2(FC) 120h+Dox. vs 48h+Dox.", y="-log10(p-value)") + 
    theme_pubr() +
    theme(strip.text.x = element_text(size=6))

plts[["dge-volcano-fc_vs_padj-t"]] = df %>% 
    ggplot(aes(x=log2FoldChange, y=(-log10_padj))) + 
    geom_scattermore(pointsize=8, alpha=0.5, pixels=c(1000,1000)) + 
    geom_point(aes(x=log2FoldChange, y=(-log10_padj)), 
               df %>% filter(gene%in%genes_oi), 
               color="orange") +
    geom_hline(yintercept=-log10(0.1), size=0.2, linetype="dashed") + 
    facet_wrap(~cell_line, scales="free_y") + 
    labs(x="log2(FC) 120h+Dox. vs 48h+Dox.", y="-log10(FDR)") + 
    theme_pubr() +
    theme(strip.text.x = element_text(size=6))

# effect of overexpressing oncogenes
df = dge_oe
plts[["dge-volcano-fc_vs_pvalue-oe"]] = df %>% 
    ggplot(aes(x=log2FoldChange, y=(-log10_pvalue))) + 
    geom_scattermore(pointsize=8, alpha=0.5, pixels=c(1000,1000)) + 
    geom_point(aes(x=log2FoldChange, y=(-log10_pvalue)), 
               df %>% filter(gene%in%genes_oi), 
               color="orange") +
    geom_hline(yintercept=-log10(0.05), size=0.2, linetype="dashed") + 
    facet_wrap(~cell_line, scales="free_y") + 
    labs(x="log2(FC) Dox. vs Ctl.", y="-log10(p-value)") + 
    theme_pubr() +
    theme(strip.text.x = element_text(size=6))

plts[["dge-volcano-fc_vs_padj-oe"]] = df %>% 
    ggplot(aes(x=log2FoldChange, y=(-log10_padj))) + 
    geom_scattermore(pointsize=8, alpha=0.5, pixels=c(1000,1000)) + 
    geom_point(aes(x=log2FoldChange, y=(-log10_padj)), 
               df %>% filter(gene%in%genes_oi), 
               color="orange") +
    geom_hline(yintercept=-log10(0.1), size=0.2, linetype="dashed") + 
    facet_wrap(~cell_line, scales="free_y") + 
    labs(x="log2(FC) Dox. vs Ctl.", y="-log10(FDR)") + 
    theme_pubr() +
    theme(strip.text.x = element_text(size=6))
    

# save
dir.create(output_figdir, recursive=TRUE)

plt_name = "metadata-eda"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 10, height = 12)

plt_name = "genexpr-eda-oe"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

plt_name = "genexpr-eda-t"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 13)

plt_name = "genexpr-eda-oe-norm"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

plt_name = "genexpr-eda-t-norm"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 13)

# TP53 mutation
plt_name = "dge-volcano-fc_vs_pvalue-mut"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

plt_name = "dge-volcano-fc_vs_padj-mut"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

# time
plt_name = "dge-volcano-fc_vs_pvalue-t"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

plt_name = "dge-volcano-fc_vs_padj-t"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

# overexpression
plt_name = "dge-volcano-fc_vs_pvalue-oe"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

plt_name = "dge-volcano-fc_vs_padj-oe"
ggsave(sprintf(file.path(output_figdir,"%s.pdf"), plt_name), plts[[plt_name]], units = 'cm', width = 12, height = 12)

print("Done!")