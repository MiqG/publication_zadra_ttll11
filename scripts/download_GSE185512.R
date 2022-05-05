# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Download files from GEO.
#

require(tidyverse)
require(GEOquery)

# variables
ROOT = here::here()
DATA_DIR = file.path(ROOT,'data')
RAW_DIR = file.path(DATA_DIR,'raw')

# inputs
GEO_ID = "GSE185512"
COUNTS_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE185512&format=file&file=GSE185512%5FExpression%5Fdataset%2Ecsv%2Egz"

# outputs
metadata_file = file.path(RAW_DIR,GEO_ID,"metadata.tsv")
counts_file = file.path(RAW_DIR,GEO_ID,"genexpr_counts.tsv.gz")

# download metadata from GEO
Sys.setenv("VROOM_CONNECTION_SIZE"=1e6)
gse = getGEO(GEO=GEO_ID)
metadata = gse[[1]] %>% phenoData() %>% pData() %>% rownames_to_column("sampleID")

# preprocess metadata
metadata = metadata %>%
    mutate(
        sample_code = gsub("\\]","",gsub("\\[","",gsub(".*? ", "", title))),
        cell_line = gsub("\\(P53-/-\\)","mutTP53",`genotype/variation:ch1`),
        treatment = `treatment:ch1`,
        timepoint = `time point:ch1`,
        replicate = gsub(".*_","",gsub("\\s+[^ ]+$", "", title))
    )

# download gene counts matrix
con = gzcon(url(COUNTS_URL))
counts = read_csv(con)

# preprocess gene counts matrix
counts = counts %>% 
    rename_at(vars(metadata[["sample_code"]]), ~ metadata[["sampleID"]]) %>%
    rename(gene = probe)

# save
dir.create(dirname(metadata_file), recursive=TRUE)
write_tsv(metadata, metadata_file)
write_tsv(counts, counts_file)

print("Done!")