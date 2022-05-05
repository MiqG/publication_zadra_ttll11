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
require(limma)
require(optparse)

# variables
SAMPLE_TYPES_OI = c('Solid Tissue Normal','Primary Tumor')

# Development
# -----------
# ROOT = here::here()
# DATA_DIR = file.path(ROOT,'data')
# RAW_DIR = file.path(DATA_DIR,'raw')
# PREP_DIR = file.path(DATA_DIR,'prep')
# phenotype_file = file.path(PREP_DIR,'sample_phenotype.tsv')
# methylation_file = file.path(PREP_DIR,"methylation_TCGA-tss.tsv.gz")
# mapping_file = file.path(PREP_DIR,"methylation_probes_tss_mapping.tsv")

# # ouptut
# output_file = file.path(PREP_DIR,'methylation_TCGA-tss-by_gene.tsv.gz')

main = function(){
    # parse arguments
    option_list = list( 
        make_option("--phenotype_file", type="character"),
        make_option("--methylation_file", type="character"),
        make_option("--mapping_file", type="character"),
        make_option("--cancer_oi", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    phenotype_file = args$phenotype_file
    methylation_file = args$methylation_file
    mapping_file = args$mapping_file
    cancer_oi = args$cancer_oi
    output_file = args$output_file
    
    # create output directory
    dir.create(dirname(output_file))
    
    # load
    print("Loading ...")
    mapping = read_tsv(mapping_file)
    metadata = read_tsv(phenotype_file) %>%
        filter(cancer_type %in% cancer_oi) %>%
        filter(sample_type %in% SAMPLE_TYPES_OI)
    
    # consider only corresponding samples
    samples_oi = metadata %>%
        filter(cancer_type %in% cancer_oi) %>%
        pull(sample)
    methylation = read.columns(methylation_file, required.col = c('sample',samples_oi))

    # drop duplicate samples
    is_duplicated = colnames(methylation)[duplicated(colnames(methylation))]
    print(sprintf("Dropping %s duplicated samples.", length(is_duplicated)))
    X = methylation %>% dplyr::select(!all_of(is_duplicated))

    # each gene may have multiple TSSs, we consider the union
    X = mapping %>%
        left_join(X, by=c("methylation_probe"="sample"))
    X = X %>%
        filter(!duplicated(methylation_probe))

    # for each gene and sample, we compute the median methylation 
    # across all probes in its united TSS
    print("Computing medians ...")
    methylation_meds = X %>%
        dplyr::select(!all_of(c("methylation_probe", "ensembl_transcript"))) %>%
        group_by(ensembl_gene, symbol) %>%
        summarize_if(is.numeric, median, na.rm=TRUE) %>%
        ungroup()
    
    # save
    print("Saving ...")
    write_tsv(methylation_meds, output_file)

}
#### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
