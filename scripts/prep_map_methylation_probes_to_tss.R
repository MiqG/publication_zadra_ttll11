# Script purpose
# --------------
# Get 
#
# Outline
# -------
# 1. Get hg19 promoter coordinates (1000 bp upstream) --> 1000TSS
# 2. 
#
# References
# ----------
#  1. https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-018-0205-1#Sec2

# load libraries
require(tidyverse)
require(GenomicRanges)

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
RAW_DIR = file.path(DATA_DIR,'raw')
PREP_DIR = file.path(DATA_DIR,'prep')

# variables
CANCERS_OI = c('LUSC','UCEC','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')

# input
tss_annotation_file = file.path(RAW_DIR,'GENCODE','gene_tss-hg19-wup1500_wdown0.tsv.gz')
probe_map_file = file.path(RAW_DIR,'UCSCXena','TCGA','methylation','probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy')
#phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
#methylation_file = file.path(RAW_DIR,'UCSCXena','TCGA','methylation','jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz')

# ouptut
output_file = file.path(PREP_DIR,'methylation_probes_tss_mapping.tsv')
probes_oi_file = file.path(PREP_DIR,'methylation_probes_tss_list.txt')

# load TSS coordinates and methylation probe map
tss = read_tsv(tss_annotation_file) %>%
    mutate(chr = paste0("chr",chr))
probes = read_tsv(probe_map_file) %>%
    dplyr::rename(chr=chrom, start=chromStart, end=chromEnd) %>%
    drop_na(`#id`,chr,start,end,strand) %>%
    mutate(strand=gsub("\\.","*",strand))

# get methylation probes in TSS1000
make_granges = function(coords, names){
    GRanges(
        seqnames = coords$chr,
        ranges = IRanges(
            start = coords$start,
            end = coords$end,
            names = coords[[names]]
        ),
        strand = coords$strand
    )
}

get_intersection = function(x,y){
    # get seqnames in x that intersect with y.
    res = findOverlaps(x,y)
    df = data.frame(
        names(x)[queryHits(res)],
        names(y)[subjectHits(res)]
    )
    colnames(df) = c('methylation_probe','ensembl_transcript')
    return(df)
}

gr_tss = make_granges(tss, "ensembl_transcript")
gr_probes = make_granges(probes, "#id")

mapping = get_intersection(gr_probes, gr_tss)

probes_tss = mapping %>%
    left_join(tss %>% distinct(ensembl_transcript, ensembl_gene, symbol), 
              by="ensembl_transcript")

probes_oi = probes_tss %>% pull(methylation_probe) %>% unique()

# save
write_tsv(probes_tss, output_file)
writeLines(probes_oi, probes_oi_file)

print('Done!')
