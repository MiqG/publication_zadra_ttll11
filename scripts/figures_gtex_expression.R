# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Visualize the expression of TTLL11 and TTLL13 across human tissues from GTEx.
# 

require(tidyverse)
require(ggpubr)
require(pheatmap)

# variables
GENES_OI = c('TTLL11','TTLL13')

ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
RESULTS_DIR = file.path(ROOT,'results')

# inputs
gene_tpms = file.path(DATA_DIR,'prep','genexpr_gtex.tsv.gz')
sampleattr = file.path(DATA_DIR,'raw','GTEx','v8','annotations','GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')

# outputs
output_fig = file.path(RESULTS_DIR,'figures','gtex-heatmap.pdf')

##### SCRIPT #####
# load data
genexpr = read_tsv(gene_tpms)
metadata = read_tsv(sampleattr)

# prep
genexpr = genexpr %>% 
    dplyr::select(-one_of('Name')) %>% 
    dplyr::rename(gene=Description) %>%
    mutate(gene=gsub('TTLL13P', 'TTLL13', gene)) %>% # exchange to popular name
    filter(gene %in% GENES_OI)

# visualize
make_heatmap = function(X, metadata, params){
    # get parameters
    genes_oi = params[['genes_oi']]
    metadata_sample_col = params[['metadata_sample_col']]
    metadata_cols_oi = params[['metadata_cols_oi']]
    palettes_cols_oi = params[['palettes_cols_oi']]
    log_transform = params[['log_transform']]
    
    # prepare matrix
    mat = X
    mat[is.na(mat)] = 0
    if ( any(genes_oi != 'NA') ){mat = mat[genes_oi,]}
    if ( log_transform ){ mat = log2(mat + 1) }
    
    # prepare annotation
    annot = metadata %>% 
        column_to_rownames(metadata_sample_col) %>%
        dplyr::select(all_of(metadata_cols_oi))
    annot = annot[(rownames(annot) %in% colnames(X)),,drop=FALSE]
    
    # prepare annotation colors
    colors = sapply(1:length(metadata_cols_oi), simplify = FALSE, 
                    function(i){
                        x = unique(annot[,metadata_cols_oi[i]])
                        palette = get_palette(palettes_cols_oi[i], length(x))
                        color = setNames(palette, x)
                        return(color)
                    })
    names(colors) = metadata_cols_oi
    
    # plot
    plts = list()
    plts[['heatmap']] = pheatmap(t(mat), 
                                 show_rownames = FALSE, 
                                 annotation_row = annot, 
                                 annotation_colors = colors, 
                                 color = get_palette('Reds',15),
                                 silent = TRUE)
    return(plts)
}

params = list(
    'genes_oi' = GENES_OI,
    'metadata_sample_col' = 'SAMPID',
    'metadata_cols_oi' = 'SMTS',
    'palettes_cols_oi' = 'Paired',
    'log_transform' = TRUE
)
plt = make_heatmap(genexpr %>% column_to_rownames('gene'), metadata, params)


# save
ggsave(output_fig, plt[['heatmap']], width=17, height=25, dpi=300, units='cm', limitsize=FALSE)


print('Done!')