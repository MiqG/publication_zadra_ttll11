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
# 1. Download data
# 2. Preprocess data:
#    - Subset phenotype data with samples of interest (>20 STN and PT samples per cancer).
#    - clean aneuploidy and centrosome amplification scores downloaded from articles.
# 3. Expression of TTLL11 and TTLL13 across human tissues.
# 4. Differential expression of TTLL11 across cancer types.
# 5. Correlation between aneuploidy score and expression of TTLL11 across cancer types.
# 6. Mutation frequency per kilobase of every gene in primary tumors.

import os
import pathlib

##### VARIABLES #####
ROOT = pathlib.Path().parent.absolute()
DATA_DIR = os.path.join(ROOT,'data')
RAW_DIR = os.path.join(DATA_DIR,'raw')
GTEX_DIR = os.path.join(RAW_DIR,'GTEx','v8')
XENA_DIR = os.path.join(RAW_DIR,'UCSCXena')
ARTICLES_DIR = os.path.join(RAW_DIR,'articles')

print(ROOT)

##### RULES ######
rule all:
    input:
        # download data
        '.done/GTEx.done', # gene expression in human tissues
        '.done/UCSCXena-TCGA-PANCAN.done', # gene expression in tumors
        '.done/Taylor2018.done', # TCGA aneuploidy

##### 1. Download data #####
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