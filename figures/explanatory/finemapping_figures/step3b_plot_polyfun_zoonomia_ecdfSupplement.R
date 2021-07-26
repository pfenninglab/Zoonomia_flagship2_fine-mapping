ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(wesanderson)

####################################
## read in table of fine-mapped SNPs ##
PROJDIR='figures/explanatory/finemapping_figures'
i_am(file.path(PROJDIR, 'step3_plot_polyfun_zoonomia_ecdf.R'))


##############################
# colors for conserved and accelerated 241mam phyloP
top_phyloP_lvls = c('Con.top.0-1%', 'Con.top.1-2%', 'Con.top.2-3%', 'Con.top.3-4%',
                    'Acc.top.3-4%', 'Acc.top.2-3%', 'Acc.top.1-2%', 'Acc.top.0-1%',  
                    'Other')
top_phyloP_cols = c(brewer.pal(4,'Blues'),brewer.pal(4,'Reds'),'#bdbdbd')
names(top_phyloP_cols) = c('Con.top.3-4%', 'Con.top.2-3%', 'Con.top.1-2%', 'Con.top.0-1%', 
                           'Acc.top.3-4%', 'Acc.top.2-3%', 'Acc.top.1-2%', 'Acc.top.0-1%',  
                           'Other')

# colors for 43primate phastCons
top_phastCons_lvls = c('PhastCons.top.0-1%', 'PhastCons.top.1-2%', 'PhastCons.top.2-3%', 
                       'PhastCons.top.3-4%', 'PhastCons.top.4-5%', 'Other')
top_phastCons_cols = c(rev(brewer.pal(5,'PuBuGn')),'#bdbdbd')
names(top_phastCons_cols) = top_phastCons_lvls

# colors for 3 ENCODE3 cCRE annotations
ENCODE3_cCRE_lvls = c('PLS', 'pELS', 'dELS', 'Other')
ENCODE3_cCRE_cols = c(brewer.pal(3,'Dark2'),'#bdbdbd')
names(ENCODE3_cCRE_cols) = ENCODE3_cCRE_lvls


########################################
## read in the fine-mapping dataframe ##
poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20210520.rds')
snps_df = readRDS(file = poly_fn)
snps_df = snps_df %>% 
  mutate(inPhyloP = case_when(grepl('Con', top_phyloP) ~ 'PhyloP_Cons', 
                              grepl('Acc', top_phyloP) ~ 'PhyloP_Acc',
                              TRUE ~ 'Not in PhyloP'),
         inPhyloP = relevel(factor(inPhyloP), ref = 'Not in PhyloP'),
         inPhastCons = case_when(grepl('Con', top_phastCons) ~ 'PhastCons', 
                                 TRUE ~ 'Not in PhastCons'),
         inPhastCons = relevel(factor(inPhastCons), ref = 'Not in PhastCons'),
         cCRE_group = relevel(cCRE_group, ref = 'Other'))

table(snps_df$group)


##########################################################
## 1. make the main ECDF plots showing PIP improvements ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_fig = 1.75; width_fig = 2.25; font_fig = 7
height_ppt = 4; width_ppt = 8;

## ecdf for phyloP annotations
pp1 = ggplot(snps_df %>% filter(grepl('non|ZoonomiaAnnot', group)),
             aes(x = PIP, color = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(~inPhyloP, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_colour_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank())  +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))

## ecdf for phastCons annotations
pp2 = ggplot(snps_df %>% filter(grepl('non|ZoonomiaAnnot', group)),
             aes(x = PIP, color = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(~inPhastCons, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_colour_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))

## ecdf for cCRE annotations
pp3 = ggplot(snps_df %>% filter(grepl('non|ZoonomiaAnnot', group)),
             aes(x = PIP, color = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(~cCRE_group, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_colour_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))

## put the eCDF panels together
plot_fn = here(PROJDIR,'plots',
               paste0('supplemental_polyfun_zoonomia_finemapping_ecdf_20210724.fig.pdf'))
pdf(plot_fn, height = 4, width = width_fig, onefile = F)
pp = ggarrange(pp1, pp2, pp3, 
               font.label = list(size = font_fig + 4, color = "black", face = "bold"),
               common.legend = TRUE, legend = 'top', ncol = 1, nrow = 3)
print(pp)
dev.off()




