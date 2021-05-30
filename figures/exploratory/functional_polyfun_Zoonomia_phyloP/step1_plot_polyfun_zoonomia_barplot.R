ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(rtracklayer)

PROJDIR='figures/exploratory/functional_polyfun_Zoonomia_phyloP'
i_am(file.path(PROJDIR, 'step1_plot_polyfun_zoonomia_barplot.R'))

poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20210520.rds')
snps_df = readRDS(file = poly_fn)

#################################
## make plots for presentation ##
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

# colors for 2 HAR and CHAR annotations
HAR_lvls = c('HAR1500bp', 'CHAR1500bp', 'Other')
HAR_cols = c(brewer.pal(3,'Accent'))
names(HAR_cols) = HAR_lvls

#################################
## make plots for presentation ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7

plot_fn = here(PROJDIR,'plots',
               paste0('polyfun_zoonomia_finemapping_20210521.ppt.pdf'))
# pdf(plot_fn, height = height_ppt, width = width_ppt)
# for (cutoff in c(.95, .9, .5, .25,0.01)){
  for (cutoff in c(0.01)){
    plot_fn = here(PROJDIR,'plots',
  paste0('polyfun_zoonomia_finemapping_PIP',cutoff,'_20210521.ppt.pdf'))
pdf(plot_fn, height = height_ppt, width = width_ppt)

pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = group)) +
  geom_bar(aes(fill = top_phyloP)) + 
  geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
  scale_fill_manual(values = top_phyloP_cols) + 
  facet_wrap(~TRAIT, scales = 'free_y', ncol = 6) + 
  scale_y_continuous(expand = expansion(mult = c(0, .4))) + 
  xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
  theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)

pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = group)) +
  geom_bar(aes(fill = top_phastCons)) + 
  geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
  scale_fill_manual(values = top_phastCons_cols) + 
  facet_wrap(~TRAIT, scales = 'free_y', ncol = 6) + 
  scale_y_continuous(expand = expansion(mult = c(0, .4))) + 
  xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
  theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)

pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = group)) +
  geom_bar(aes(fill = cCRE_group)) + 
  geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
  scale_fill_manual(values = ENCODE3_cCRE_cols) + 
  facet_wrap(~TRAIT, scales = 'free_y', ncol = 6) + 
  scale_y_continuous(expand = expansion(mult = c(0, .4))) + 
  xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
  theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)

dev.off()
}
