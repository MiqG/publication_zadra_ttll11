# Differential gene expression analysis
# 1. load gene counts, metadata, design formula
# 2. filter out low expressed genes; keep genes with at least 10 reads total
# 3. do DESeq
# 4. save result

require(optparse)
require(tidyverse)
require(DESeq2)

# variables
MIN_GENE_READS = 10

# Development
# -----------
# ROOT = here::here()
# DATA_DIR = file.path(ROOT,'data')
# RAW_DIR = file.path(DATA_DIR,'raw')
# GEO_ID = "GSE185512"
# metadata_file = file.path(RAW_DIR,GEO_ID,"metadata.tsv")
# counts_file = file.path(RAW_DIR,GEO_ID,"genexpr_counts.tsv.gz")
# sample_col = "sampleID"
# subset_cols = "cell_line|timepoint|treatment"
# subset_values = "RPE-1-Empty;RPE-1 (P53-/-)-Empty|120h|NO Dox;Doxycycline"
# design = "~ timepoint + treatment + cell_line"
# comparison_col = "cell_line"
# condition_a = "RPE-1 (P53-/-)-Empty"
# condition_b = "RPE-1-Empty"


##### FUNCTIONS #####
load_data = function(counts_file, metadata_file, sample_col, subset_cols, subset_values){
    # gene expression
    counts = as.data.frame(read_tsv(counts_file))
    # create index
    rownames(counts) = counts[[1]]
    counts = counts[,-1]
    
    # metadata
    metadata = as.data.frame(read_tsv(metadata_file)) %>% 
        column_to_rownames(sample_col)
    
    # filter samples of interest
    subset_cols_split = strsplit(subset_cols, "\\|")[[1]]
    subset_values_split = strsplit(subset_values, "\\|")[[1]]
    subsets = setNames(subset_values_split, subset_cols_split)
    for (subset_col in names(subsets)){
        subset_values_vec = strsplit(subsets[[subset_col]],";")[[1]]
        metadata = metadata %>% 
            filter(metadata[[subset_col]] %in% subset_values_vec)
    }
    
    # ensure match
    common_samples = intersect(rownames(metadata),colnames(counts))
    counts = counts[,common_samples]
    metadata = metadata[common_samples,]
    
    dat = list(counts=counts, metadata=metadata)
    
    return(dat)
}


filter_genes = function(counts){
    rsums = rowSums(counts)
    counts = counts[rsums >= MIN_GENE_READS,]
    return(counts)
}


create_dataset = function(counts, metadata, design){
    ds = DESeqDataSetFromMatrix(counts, metadata, design)
    return(ds)
}


run_DESeq = function(counts, metadata, design, comparison_col, condition_a, condition_b){
    result = tryCatch({
        # prepare inputs
        ds = create_dataset(counts, metadata, design)
        # run
        ds = DESeq(ds)
        # prepare outputs
        result = results(ds, contrast=c(comparison_col, condition_a, condition_b))
        result = cbind(gene=result@rownames, as.data.frame(result))
        result = result %>%
            mutate(
                log10_pvalue = log10(pvalue),
                log10_padj = log10(padj),
                comparison_col = comparison_col,
                condition_a = condition_a,
                condition_b = condition_b
            )
        return(result)
    },
      error = function(cond){
          result = data.frame(
              gene=rownames(counts),
              baseMean=NA,
              log2FoldChange=NA,
              lfcSE=NA,
              stat=NA,
              pvalue=NA,
              padj=NA,
              log10_pvalue=NA,
              log10_padj=NA
          )
          result = result %>%
            mutate(
                log10_pvalue = log10(pvalue),
                log10_padj = log10(padj),
                comparison_col = comparison_col,
                condition_a = condition_a,
                condition_b = condition_b
            )
          print(paste('[Catched eror]:',cond))
          return(result)
      })
    
    
    return(result)
}


main = function(){
    # parse arguments
    option_list = list( 
        make_option("--counts_file", type="character"),
        make_option("--metadata_file", type="character"),
        make_option("--design", type="character"),
        make_option("--sample_col", type="character"),
        make_option("--subset_cols", type="character"),
        make_option("--subset_values", type="character"),
        make_option("--comparison_col", type="character"),
        make_option("--condition_a", type="character"),
        make_option("--condition_b", type="character"),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    counts_file = args$counts_file
    metadata_file = args$metadata_file
    design = args$design
    sample_col = args$sample_col
    subset_cols = args$subset_cols
    subset_values = args$subset_values
    comparison_col = args$comparison_col
    condition_a = args$condition_a
    condition_b = args$condition_b
    output_file = args$output_file
    
    # create output directory
    dir.create(dirname(output_file))
    
    # load
    dat = load_data(counts_file, metadata_file, sample_col, 
                    subset_cols, subset_values)
    counts = filter_genes(dat[["counts"]])
    metadata = dat[["metadata"]]
    
    # run
    result = run_DESeq(counts, metadata, as.formula(design), 
                       comparison_col, condition_a, condition_b)
    result = result %>%
        mutate(design = design,
               sample_col = sample_col,
               subset_cols = subset_cols,
               subset_values = subset_values)
    print(dim(result))
    
    # save
    print(sprintf("Saving at %s ...", output_file))
    write_tsv(result, output_file)
}


# #### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
