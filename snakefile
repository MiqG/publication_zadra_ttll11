# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Workflow to carry out whole bioinformatic analysis for the article 
# Zadra *et al.* 2021 (DOI: ).
#
# Outline
# -------
# 1. Download of publicly available datasets from TCGA:
#      - PANCAN-normalized gene expression
#      - PANCAN Copy Number Variation (CNV)
#      - PANCAN Somatic Nucleotide Variation (SNV)
#      - PANCAN Aneuploidy Scores
#      - PANCAN Centrosome Amplification
#      - PANCAN Phenotype Data
# 2. Preprocess data:
#    - Subset phenotype data with samples of interest (>20 STN and PT samples per cancer).
#    - clean aneuploidy and centrosome amplification scores downloaded from articles.
# 3. Assess and plot differential expression of *TTLL11* by cancer type.
# 4. Plot changes of CNV by cancer type.
# 5. Measure mutation type frequency per gene per kilobase per sample across cancer type.
# 6. Plot measured mutation frequencies by cancer type.
# 7. Measure spearman correlation between aneuploidy scores and expression of all genes by cancer type.
# 9. Measure spearman correlation between centrosome amplification scores and expression of all genes by cancer type.
# 10. Plot measured correlations with aneuploidy and centrosome scores. 

##### VARIABLES #####

##### RULES ######
