# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Workflow purpose
# ----------------
# Carry out whole bioinformatic analysis for the article:
# Zadra *et al.* XXXX (DOI: XXXX).
#
# Outline
# -------
# 0. Download data
# 1. Preprocess data:
#    - Subset phenotype data with samples of interest (>20 STN and PT samples per cancer).
#    - clean aneuploidy scores downloaded from article.
# 2. Expression of TTLL11 and TTLL13 across human tissues.
# 3. Differential expression of TTLL11 across cancer types.
# 4. Correlation between aneuploidy score and expression of TTLL11 across cancer types.
# 5. Mutation frequency per kilobase of every gene in primary tumors.
# 6. Prepare publishable figures

import os
import pathlib

##### VARIABLES #####
ROOT = pathlib.Path().parent.absolute()
DATA_DIR = os.path.join(ROOT,'data')
RAW_DIR = os.path.join(DATA_DIR,'raw')
PREP_DIR = os.path.join(DATA_DIR,'prep')
GTEX_DIR = os.path.join(RAW_DIR,'GTEx','v8')
XENA_DIR = os.path.join(RAW_DIR,'UCSCXena')
ARTICLES_DIR = os.path.join(RAW_DIR,'articles')
RESULTS_DIR = os.path.join(ROOT,'results')

##### RULES ######
rule all:
    input:
        # download data
        '.done/GTEx.done', # gene expression in human tissues
        '.done/UCSCXena-TCGA-PANCAN.done', # gene expression in tumors
        '.done/Taylor2018.done', # TCGA aneuploidy
        
        # preprocess data
        '.done/prep-sample_phenotype.done',
        '.done/prep-aneuploidy.done',
        '.done/prep-snv.done',
        '.done/prep-genexpr_TTLL11.done',
        '.done/prep-gtex.done',
        
        # expression of TTLL11 and TTLL13 in tissues
        '.done/genexpr_human_tissues.done',
        
        # differential expression of TTLL11 across cancer types
        '.done/diffexpr_TTLL11.done',
        
        # correlation gene expression vs aneuploidy score
        '.done/aneuploidy_correlation.done',
        
        # mutation frequency
        '.done/mutation_frequency.done',
        
        # publish figures
        '.done/publish_figures.done'
        
        

##### 0. Download data #####
rule download_gtex:
    params:
        # annotations
        sampleattr_desc = 'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx',
        subjpheno_desc = 'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx',
        sampleattr = 'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
        subjpheno = 'https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt',
        
        # RNA-Seq
        gene_tpms = 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz'
        
    output:
        touch('.done/GTEx.done'),
        
        # annotations
        sampleattr_desc = os.path.join(GTEX_DIR,'annotations','GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx'),
        subjpheno_desc = os.path.join(GTEX_DIR,'annotations','GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx'),
        sampleattr = os.path.join(GTEX_DIR,'annotations','GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'),
        subjpheno = os.path.join(GTEX_DIR,'annotations','GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'),
        
        # RNA-Seq
        gene_tpms = os.path.join(GTEX_DIR,'rnaseq','GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
        
    shell:
        """
        # annotations
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.sampleattr_desc} -O {output.sampleattr_desc}
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.subjpheno_desc} -O {output.subjpheno_desc}
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.sampleattr} -O {output.sampleattr}
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.subjpheno} -O {output.subjpheno}
        
        # RNA-Seq
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.gene_tpms} -O {output.gene_tpms}

        echo Done!
        """
        
        
rule download_ucscxena_tcga_pancan:
    params:
        # RNA-Seq
        genexpr = 'https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz',
        
        # Mutations
        snv = 'https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/mc3.v0.2.8.PUBLIC.xena.gz',
        
        
        # Phenotype
        clinical = 'https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp',
        sample_type ='https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz'
        
    output:
        touch('.done/UCSCXena-TCGA-PANCAN.done'),
        
        # RNA-Seq
        genexpr = os.path.join(XENA_DIR,'TCGA','rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz'),
        
        # Mutations
        snv = os.path.join(XENA_DIR,'TCGA','snv','mc3.v0.2.8.PUBLIC.xena.gz'),
        
        # phenotype
        clinical = os.path.join(XENA_DIR,'TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz'),
        sample_type = os.path.join(XENA_DIR,'TCGA','phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz'),
        
        readme = os.path.join(XENA_DIR,'TCGA','README.md')
    
    shell:
        """
        # RNA-Seq
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.genexpr} -O {output.genexpr}
        
        # Mutations
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.snv} -O {output.snv}
        
        # Phenotype
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.clinical} -O {output.clinical}
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.sample_type} -O {output.sample_type}        
        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """

rule download_tcga_aneuploidy_scores:
    message:
        "Aneuploidy scores from Taylor *et al.* 2018 (https://doi.org/10.1016/j.ccell.2018.03.007)."
    params:
        aneuploidy = 'https://www.cell.com/cms/10.1016/j.ccell.2018.03.007/attachment/2d887978-0a9c-4e90-af00-66eb5948673b/mmc2.xlsx'
    output:
        touch('.done/Taylor2018.done'),
        aneuploidy = os.path.join(ARTICLES_DIR,'Taylor2018','aneuploidy.xlsx'),
        readme = os.path.join(ARTICLES_DIR,'Taylor2018','README.md')
    shell:
        """
        # aneuploidy
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.aneuploidy} -O {output.aneuploidy}        
        
        # readme
        echo "Downloaded on $(date)." > {output.readme}
        
        echo Done!
        """ 
        
##### 1. Preprocess data #####
rule prep_sample_phenotype:
    input:
        clinical = os.path.join(XENA_DIR,'TCGA','phenotype','Survival_SupplementalTable_S1_20171025_xena_sp.gz'),
        sample_type = os.path.join(XENA_DIR,'TCGA','phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz')
    output:
        touch('.done/prep-sample_phenotype.done'),
        os.path.join(PREP_DIR,'sample_phenotype.tsv')
    shell:
        """
        Rscript scripts/prep_sample_phenotype.R
        """

rule prep_aneuploidy:
    input:
        aneuploidy = os.path.join(ARTICLES_DIR,'Taylor2018','aneuploidy.xlsx')
    output:
        touch('.done/prep-aneuploidy.done'),
        os.path.join(PREP_DIR,'aneuploidy.tsv')
    shell:
        """
        Rscript scripts/prep_aneuploidy_scores.R
        """
        
        
rule prep_snv:
    input:
        snv = os.path.join(XENA_DIR,'TCGA','snv','mc3.v0.2.8.PUBLIC.xena.gz'),
        sample_phenotype = os.path.join(PREP_DIR,'sample_phenotype.tsv')
    output:
        touch('.done/prep-snv.done'),
        os.path.join(PREP_DIR,'snv_gene_freq.tsv')
    shell:
        """
        Rscript scripts/prep_snv_counts.R
        """

        
rule prep_genexpr_TTLL11:
    input:
        genexpr = os.path.join(XENA_DIR,'TCGA','rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz'),
        sample_phenotype = os.path.join(PREP_DIR,'sample_phenotype.tsv')
    output:
        touch('.done/prep-genexpr_TTLL11.done'),
        os.path.join(PREP_DIR,'genexpr_TTLL11.tsv')
    shell:
        """
        Rscript scripts/prep_TTLL11_expression.R
        """

rule prep_gtex:
    input:
        gene_tpms = os.path.join(GTEX_DIR,'rnaseq','GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
    output:
        touch('.done/prep-gtex.done'),
        gene_tpms = os.path.join(PREP_DIR,'genexpr_gtex.tsv.gz')
    shell:
        """
        zgrep -e "GTEX\|TTLL11\|TTLL13" {input} | gzip > {output.gene_tpms}
        """
        
    
##### 2. Gene expression of TTLL11 and TTLL13 across human tissues #####
rule gtex_genexpr:
    input:
        gene_tpms = os.path.join(PREP_DIR,'genexpr_gtex.tsv.gz'),
        sampleattr = os.path.join(GTEX_DIR,'annotations','GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt') 
    output:
        touch('.done/genexpr_human_tissues.done'),
        os.path.join(RESULTS_DIR,'figures','gtex-heatmap.pdf')
    shell:
        """
        Rscript scripts/figures_gtex_expression.R
        """
        
##### 3. Differential expression of TTLL11 across cancer types #####
rule diffexpr_TTLL11:
    input:
        genexpr = os.path.join(PREP_DIR,'genexpr_TTLL11.tsv'),
    output:
        touch('.done/diffexpr_TTLL11.done'),
        rds = os.path.join(RESULTS_DIR,'figures','differential_expression.rds'),
        figdata = os.path.join(RESULTS_DIR,'files','figdata-differential_expression.xlsx')
    shell:
        """
        Rscript scripts/figures_differential_genexpr.R
        """

##### 4. Correlation between aneuploidy score and expression of TTLL11 #####
rule correlate_aneuploidy:
    input:
        aneuploidy_scores = os.path.join(PREP_DIR,'aneuploidy.tsv'),
        phenotype = os.path.join(PREP_DIR,'sample_phenotype.tsv'),
        genexpr = os.path.join(XENA_DIR,'TCGA','rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')
    threads: 10
    output:
        correlations = os.path.join(RESULTS_DIR,'files','correlation-genexpr_aneuploidy.tsv')
    shell:
        """
        Rscript scripts/cor_genexpr_aneuploidy.R
        """
        
rule figures_gene_aneuploidy_correlation:
    input:
        correlations = os.path.join(RESULTS_DIR,'files','correlation-genexpr_aneuploidy.tsv')
    output:
        touch('.done/aneuploidy_correlation.done'),
        rds = os.path.join(RESULTS_DIR,'figures','correlation_with_scores.rds'),
        figdata = os.path.join(RESULTS_DIR,'files','figdata-correlation_with_scores.xlsx')
    shell:
        """
        Rscript scripts/figures_aneuploidy_correlation.R
        """

##### 5. Mutation frequency per kilobase of every gene in primary tumors #####
rule mutation_frequency:
    input:
        snv_freq = os.path.join(PREP_DIR,'snv_gene_freq.tsv')
    output:
        touch('.done/mutation_frequency.done'),
        rds = os.path.join(RESULTS_DIR,'figures','mutation_frequency.rds'),
        figdata = os.path.join(RESULTS_DIR,'files','figdata-mutation_frequency.xlsx')
    shell:
        """
        Rscript scripts/figures_mutation_frequencies.R
        """

##### 6. Prepare figures for publication #####
rule publish_figures:
    input:
        os.path.join(RESULTS_DIR,'figures','mutation_frequency.rds'),
        os.path.join(RESULTS_DIR,'figures','differential_expression.rds'),
        os.path.join(RESULTS_DIR,'figures','correlation_with_scores.rds')
    output:
        touch('.done/publish_figures.done'),
        os.path.join(RESULTS_DIR,'figures','expression_aneuploidy.pdf'),
        os.path.join(RESULTS_DIR,'figures','mutations.pdf')
    shell:
        """
        Rscript scripts/publish_figures.R
        """
        
