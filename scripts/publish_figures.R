# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Combine plots into figures ready to publish.
#
# Outline
# -------
# - figure with TTLL11 gene expression and correlation with aneuploidy
# - figure with TTLL11 mutational analysis


require(ggpubr)
require(latex2exp)
require(patchwork)


# variables
ROOT = here::here()
RESULTS_DIR = file.path(ROOT,'results')
FIGURES_DIR = file.path(RESULTS_DIR,'figures')

# figure specs
# - variables in italics
# - type on top of color bold face
# - Capital letters for multipar figures
WIDTH_1COL = 5.5 # cm
WIDTH_2COL = 12 # cm
FONT_SIZE = 7 # pt
FONT_FAMILY = 'helvetica'
FONT_FACE = 'plain'
FIG_LETTER_FONT_SIZE = 9 # pt
FIG_LETTER_FONT_FACE = 'bold'

# inputs
filenames = c('differential_expression','correlation_with_scores','mutation_frequency')
rds_files = file.path(FIGURES_DIR,paste0(filenames,'.rds'))

# outputs
fig_expression_aneuploidy_file = file.path(FIGURES_DIR,'expression_aneuploidy')
fig_mutations_file = file.path(FIGURES_DIR,'mutations.png')

##### SCRIPT #####
# load data
plts = lapply(rds_files, function(x){ readRDS(x) })
names(plts) = filenames

##### TTLL11 expression and correlation with aneuploidy #####
# compose figure
widths_cancer_type = c(1.15,-0.1,0.13)

expression = ggarrange(
    plts[['differential_expression']][['bycancer']] + 
        ylab(TeX('$log_2(Norm. Count + 1)$')) + 
        theme_pubr(base_size = FONT_SIZE),
    NULL,
    plts[['differential_expression']][['pancancer']] + 
        theme_pubr(base_size = FONT_SIZE),
    widths = widths_cancer_type,
    common.legend = TRUE,
    ncol = 3
) %>% annotate_figure(bottom = text_grob('Cancer Type', 
                                         size = FONT_SIZE, 
                                         family = FONT_FAMILY,
                                         vjust = -0.75)
                     )

aneuploidy = ggarrange(
    plts[['correlation_with_scores']][['aneuploidy_bycancer']] + 
        ylab('Spearman Correlation Coeff.') + 
        theme_pubr(base_size = FONT_SIZE),
    NULL,
    plts[['correlation_with_scores']][['aneuploidy_pancancer']] + 
        theme_pubr(base_size = FONT_SIZE),
    widths = widths_cancer_type,
    common.legend = TRUE,
    ncol = 3
) %>% annotate_figure(bottom = text_grob('Cancer Type', 
                                         size = FONT_SIZE, 
                                         family = FONT_FAMILY,
                                         vjust = -0.75)
                     ) 

fig = expression / aneuploidy + 
    plot_annotation(tag_levels = 'A') & 
    theme(
        plot.tag = element_text(
            size=FIG_LETTER_FONT_SIZE,
            face=FIG_LETTER_FONT_FACE,
            family=FONT_FAMILY)
    )

# save
ggsave(paste0(fig_expression_aneuploidy_file,'.png'), fig, width = WIDTH_2COL, units = 'cm', dpi = 300)
ggsave(paste0(fig_expression_aneuploidy_file,'.pdf'), fig, width = WIDTH_2COL, units = 'cm', dpi = 300, device = cairo_pdf)

##### TTLL11 mutations #####
plt_names_oi = c('Missense_Mutation','Nonsense_Mutation','3\'UTR','Splice_Site','Frame_Shift_Del')
plt_names_oi = paste0('effect_by_cancer-',plt_names_oi)
mutations = sapply(plt_names_oi, function(plt_name){
    plt = plts[['mutation_frequency']][[plt_name]]
    plt = ggarrange(
        plt[['bycancer']] + 
            labs(y = TeX('Normalized $log_{10}$(Mut. Freq.)')) +
            ggtitle(gsub('effect_by_cancer-','',plt_name)) +
            theme_pubr(base_size = FONT_SIZE),
        NULL,
        plt[['pancancer']] + 
            theme_pubr(base_size = FONT_SIZE),
        widths = widths_cancer_type,
        common.legend = TRUE,
        ncol = 3
    ) %>% annotate_figure(bottom = text_grob('Cancer Type', 
                                             size = FONT_SIZE, 
                                             family = FONT_FAMILY,
                                             vjust = -0.75)
                          )
    
    return(plt)
}, simplify=FALSE)


fig = wrap_plots(mutations, ncol=1) + 
    plot_annotation(tag_levels = 'A') & 
    theme(
        plot.tag = element_text(
            size=FIG_LETTER_FONT_SIZE,
            face=FIG_LETTER_FONT_FACE,
            family=FONT_FAMILY)
    )

# save
ggsave(fig_mutations_file, fig, width = WIDTH_2COL, 
       height = 4*length(plts[['mutation_frequency']]), units = 'cm', dpi = 300)

print('Done!')