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

CANCERS_OI = [
    'LUSC',
    'UCEC',
    'BRCA',
    'STAD',
    'LUAD',
    'KIRP',
    'THCA',
    'KICH',
    'COAD',
    'LIHC',
    'HNSC',
    'PRAD',
    'KIRC'
]

CELL_LINES = [
    "RPE-1-Cdc25a",
    "RPE-1-CyclinE",
    "RPE-1-Myc",
    "RPE-1-mutTP53-Cdc25a",
    "RPE-1-mutTP53-CyclinE",
    "RPE-1-mutTP53-Myc",
    "BT549-Cdc25a",
    "BT549-CyclinE",
    "BT549-Myc",
    "HCC1806-Cdc25a",
    "HCC1806-CyclinE",
    "HCC1806-Myc",
    "MDA-MB231-Cdc25a",
    "MDA-MB231-CyclinE",
    "MDA-MB231-Myc"
]

COMPARISONS = ["dox_vs_ctl","dox_120h_vs_48h"]

##### RULES ######
rule all:
    input:
        # download data
        '.done/GTEx.done', # gene expression in human tissues
        '.done/UCSCXena-TCGA-PANCAN.done', # gene expression in tumors
        '.done/Taylor2018.done', # TCGA aneuploidy
        '.done/GSE185512.done',
        '.done/CHEA_TFs.done',
        
        # preprocess data
        '.done/prep-sample_phenotype.done',
        '.done/prep-aneuploidy.done',
        '.done/prep-snv.done',
        '.done/prep-genexpr_TTLLs.done',
        '.done/prep-gtex.done',
        os.path.join(PREP_DIR,'methylation_probes_tss_mapping.tsv'),
        os.path.join(PREP_DIR,'methylation_probes_tss_list.txt'),
        os.path.join(PREP_DIR,"methylation_TCGA-tss.tsv.gz"),
        expand(os.path.join(PREP_DIR,"methylation_TCGA-tss-by_gene",'{cancer}.tsv.gz'), cancer=CANCERS_OI),
        os.path.join(PREP_DIR,"methylation_TCGA-tss-by_gene",'merged.tsv.gz'),
        
        # expression of TTLL11 and TTLL13 in tissues
        '.done/genexpr_human_tissues.done',
        
        # differential expression of TTLL11 across cancer types
        os.path.join(RESULTS_DIR,'figures','differential_expression.pdf'),
        
        # correlation gene expression vs aneuploidy score
        os.path.join(RESULTS_DIR,'figures','aneuploidy_correlations.pdf'),
        
        # mutation frequency
        os.path.join(RESULTS_DIR,'figures','mutation_frequency.pdf'),
        
        # differential analyses
        ## run
        expand(os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{cell_line}-dox_vs_ctl.tsv"), cell_line=CELL_LINES),
        expand(os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{cell_line}-dox_120h_vs_48h.tsv"), cell_line=CELL_LINES),
        os.path.join(RESULTS_DIR,"files","dge_p53mutation-GSE185512-rpe1_120h_mut_vs_wt.tsv"),
        ## combine
        expand(os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{comparison}-merged.tsv"), comparison=COMPARISONS)
        ## visualize
        

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
        
        
rule download_GSE185512:
    output:
        touch('.done/GSE185512.done'),
        metadata = os.path.join(RAW_DIR,"GSE185512","metadata.tsv"),
        counts = os.path.join(RAW_DIR,"GSE185512","genexpr_counts.tsv.gz")
    shell:
        """
        Rscript scripts/download_GSE185512.R
        """
        
        
rule download_harmonizome_gene_sets:
    params:
        chea = "https://maayanlab.cloud/static/hdfs/harmonizome/data/cheappi/gene_set_library_crisp.gmt.gz"
    output:
        touch('.done/CHEA_TFs.done'),
        chea = os.path.join(RAW_DIR,"Harmonizome","CHEA-TranscriptionFactorTargets.gmt.gz"),
        readme = os.path.join(RAW_DIR,'Harmonizome','README.md')
    shell:
        """
        # aneuploidy
        wget --user-agent="Chrome" --no-clobber --no-check-certificate {params.chea} -O {output.chea}        
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

        
rule prep_genexpr_TTLLs:
    input:
        genexpr = os.path.join(XENA_DIR,'TCGA','rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz'),
        sample_phenotype = os.path.join(PREP_DIR,'sample_phenotype.tsv')
    output:
        touch('.done/prep-genexpr_TTLLs.done'),
        os.path.join(PREP_DIR,'genexpr_TTLLs.tsv')
    shell:
        """
        Rscript scripts/prep_TTLLs_expression.R
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
        

rule prep_methylation_tss:
    input:
        tss = os.path.join(RAW_DIR,'GENCODE','gene_tss-hg19-wup1500_wdown0.tsv.gz'),
        probe_map = os.path.join(RAW_DIR,'UCSCXena','TCGA','methylation','probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy')
    output:
        tss_gene_mapping = os.path.join(PREP_DIR,'methylation_probes_tss_mapping.tsv'),
        tss_probes_list = os.path.join(PREP_DIR,'methylation_probes_tss_list.txt')
    threads: 1
    resources:
        runtime = 3600*3, # 6h
        memory = 20 # G
    shell:
        """
        Rscript scripts/prep_map_methylation_probes_to_tss.R        
        """
        
rule tcga_subset_methylation:
    input:
        methylation = os.path.join(XENA_DIR,'TCGA','methylation','jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz'),
        probes_oi = os.path.join(PREP_DIR,'methylation_probes_tss_list.txt')
    output:
        methylation = os.path.join(PREP_DIR,"methylation_TCGA-tss.tsv.gz")
    threads: 1
    resources:
        runtime = 3600*3, # 6h
        memory = 10 # G
    shell:
        """
        set -eo pipefail
        
        zgrep --regexp="sample" --file={input.probes_oi} {input.methylation} | gzip > {output.methylation}
        
        echo "Done!"
        """

rule tcga_summarize_methylation:
    input:
        methylation = os.path.join(PREP_DIR,"methylation_TCGA-tss.tsv.gz"),
        phenotype = os.path.join(PREP_DIR,'sample_phenotype.tsv'),
        mapping = os.path.join(PREP_DIR,"methylation_probes_tss_mapping.tsv")
    output:
        os.path.join(PREP_DIR,"methylation_TCGA-tss-by_gene",'{cancer}.tsv.gz')
    params:
        cancer_oi = "{cancer}"
    threads: 1
    resources:
        runtime = 3600*3, # 6h
        memory = 20 # G
    shell:
        """
        Rscript scripts/prep_methylation_tss_tcga.R \
                    --methylation_file={input.methylation} \
                    --phenotype_file={input.phenotype} \
                    --mapping_file={input.mapping} \
                    --cancer_oi={params.cancer_oi} \
                    --output_file={output}
        """
    
rule tcga_combine_methylation:
    input:
        methylation = [os.path.join(PREP_DIR,"methylation_TCGA-tss-by_gene",'{cancer}.tsv.gz').format(cancer=cancer) for cancer in CANCERS_OI]
    output:
        methylation = os.path.join(PREP_DIR,"methylation_TCGA-tss-by_gene",'merged.tsv.gz')
    threads: 1
    resources:
        runtime = 3600*1, # 1h
        memory = 20 # G
    run:
        import os
        import pandas as pd
        
        dfs = []
        for f in input.methylation:
            print("Loading %s..." % f)
            
            df = pd.read_table(f, index_col=[0,1])
            dfs.append(df)
        
        dfs = pd.concat(dfs, axis=1)
        print("Saving...")
        dfs.reset_index().to_csv(output.methylation, sep="\t", compression="gzip", index=None)
        
        print("Done!")
            
    
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
        
##### 3. Differential expression of TTLLs across cancer types #####
rule diffexpr_TTLLs:
    input:
        genexpr = os.path.join(PREP_DIR,'genexpr_TTLLs.tsv'),
    output:
        figures = [os.path.join(RESULTS_DIR,'figures','differential_expression.pdf'), 
                   os.path.join(RESULTS_DIR,'figures','TTLLs-differential_gene_expression-barplot.pdf')],
        figdata = os.path.join(RESULTS_DIR,'files','figdata-differential_expression.xlsx')
    shell:
        """
        Rscript scripts/figures_differential_genexpr_TTLLs.R
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
        figure = os.path.join(RESULTS_DIR,'figures','aneuploidy_correlations.pdf'),
        figdata = os.path.join(RESULTS_DIR,'files','figdata-aneuploidy_correlations.xlsx')
    shell:
        """
        Rscript scripts/figures_aneuploidy_correlation.R
        """

##### 5. Mutation frequency per kilobase of every gene in primary tumors #####
rule mutation_frequency:
    input:
        snv_freq = os.path.join(PREP_DIR,'snv_gene_freq.tsv')
    output:
        figure = os.path.join(RESULTS_DIR,'figures','mutation_frequency.pdf'),
        figdata = os.path.join(RESULTS_DIR,'files','figdata-mutation_frequency.xlsx')
    shell:
        """
        Rscript scripts/figures_mutation_frequencies.R
        """

##### 6. Differential expresison in cell lines over expressing TTLL11 #####
rule run_dge_overexpression_dox_vs_ctl:
    input:
        metadata = os.path.join(RAW_DIR,"GSE185512","metadata.tsv"),
        counts = os.path.join(RAW_DIR,"GSE185512","genexpr_counts.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{cell_line}-dox_vs_ctl.tsv")
    params:
        sample_col = "sampleID",
        subset_cols = "cell_line|treatment",
        subset_values = "{cell_line}|NO Dox;Doxycycline",
        design = "~ timepoint + treatment",
        comparison_col = "treatment",
        condition_a = "Doxycycline",
        condition_b = "NO Dox"
    shell:
        """
        nice Rscript scripts/run_DESeq.R \
                    --metadata_file={input.metadata} \
                    --counts_file={input.counts} \
                    --output_file='{output}' \
                    --sample_col='{params.sample_col}' \
                    --subset_cols='{params.subset_cols}' \
                    --subset_values='{params.subset_values}' \
                    --design='{params.design}' \
                    --comparison_col='{params.comparison_col}' \
                    --condition_a='{params.condition_a}' \
                    --condition_b='{params.condition_b}'
        """
        
        
rule run_dge_overexpression_dox_120h_vs_48h:
    input:
        metadata = os.path.join(RAW_DIR,"GSE185512","metadata.tsv"),
        counts = os.path.join(RAW_DIR,"GSE185512","genexpr_counts.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{cell_line}-dox_120h_vs_48h.tsv")
    params:
        sample_col = "sampleID",
        subset_cols = "cell_line|treatment",
        subset_values = "{cell_line}|Doxycycline",
        design = "~ timepoint",
        comparison_col = "timepoint",
        condition_a = "120h",
        condition_b = "48h"
    shell:
        """
        nice Rscript scripts/run_DESeq.R \
                    --metadata_file={input.metadata} \
                    --counts_file={input.counts} \
                    --output_file='{output}' \
                    --sample_col='{params.sample_col}' \
                    --subset_cols='{params.subset_cols}' \
                    --subset_values='{params.subset_values}' \
                    --design='{params.design}' \
                    --comparison_col='{params.comparison_col}' \
                    --condition_a='{params.condition_a}' \
                    --condition_b='{params.condition_b}'
        """
        
        
rule combine_dge:
    input:
        dges = [os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{cell_line}-{comparison}.tsv").format(cell_line=cell_line, comparison="{comparison}") for cell_line in CELL_LINES]
    output:
        dges = os.path.join(RESULTS_DIR,"files","dge_overexpression-GSE185512-{comparison}-merged.tsv")
    run:
        import pandas as pd
        
        merged = pd.concat([pd.read_table(f) for f in input.dges])
        merged.to_csv(output.dges, sep="\t", index=None)
        
        print("Done!")
        
        
rule run_dge_rpe1mutation:
    input:
        metadata = os.path.join(RAW_DIR,"GSE185512","metadata.tsv"),
        counts = os.path.join(RAW_DIR,"GSE185512","genexpr_counts.tsv.gz")
    output:
        os.path.join(RESULTS_DIR,"files","dge_p53mutation-GSE185512-rpe1_120h_mut_vs_wt.tsv")
    params:
        sample_col = "sampleID",
        subset_cols = "cell_line|treatment|timepoint",
        subset_values = "RPE-1-Empty;RPE-1 mutTP53-Empty|NO Dox;Doxycycline|120h",
        design = "~ treatment + cell_line",
        comparison_col = "cell_line",
        condition_a = "RPE-1 mutTP53-Empty",
        condition_b = "RPE-1-Empty"
    shell:
        """
        nice Rscript scripts/run_DESeq.R \
                    --metadata_file={input.metadata} \
                    --counts_file={input.counts} \
                    --output_file='{output}' \
                    --sample_col='{params.sample_col}' \
                    --subset_cols='{params.subset_cols}' \
                    --subset_values='{params.subset_values}' \
                    --design='{params.design}' \
                    --comparison_col='{params.comparison_col}' \
                    --condition_a='{params.condition_a}' \
                    --condition_b='{params.condition_b}'
        
        """
        
##### 7. Methylation in TCGA #####
# rule figures_methylation_tcga:
#     input:
#         phenotype = os.path.join(),
#         methylation = os.path.join(),
#         genexpr = os.path.join()
#     output:
        
#     shell:
#         """
#         Rscript scripts/figures_methylation_tcga.R
#         """