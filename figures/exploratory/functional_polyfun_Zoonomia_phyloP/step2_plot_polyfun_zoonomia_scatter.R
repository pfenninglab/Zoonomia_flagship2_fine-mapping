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

####################################
## read in table of fine-mapped SNPs ##
PROJDIR='figures/exploratory/functional_polyfun_Zoonomia_phyloP'
i_am(file.path(PROJDIR, 'step2_plot_polyfun_zoonomia_scatter.R'))

poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20210513.rds')
snps_df = readRDS(file = poly_fn)

snps_df2 = snps_df %>%
  group_by(TRAIT, SNP) %>% filter(n() == 4,any('baselineLF2.2.UKB' %in% group)) %>%
  mutate(PIP_baseline = PIP[which(group == 'baselineLF2.2.UKB')]) %>% 
  ungroup() %>% filter(PIP_baseline > 0.01) %>%
  filter(! group %in% c( 'baselineLF2.2.UKB'))

#################################
## make plots for presentation ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7


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

plot_fn = here(PROJDIR,'plots',
               paste0('polyfun_zoonomia_finemapping_PIPscatter_20210521.ppt.pdf'))
pdf(plot_fn, height = height_ppt, width = width_ppt)
pp = ggplot(data = snps_df2, 
            aes(x = PIP_baseline, y = PIP)) + 
  # geom_point(pch = 20, alpha = .25) +
  geom_bin2d(bins = 50) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  # scale_colour_manual(values = phyloP_col) + 
  # guides(colour = guide_legend(nrow = 1)) + 
  facet_grid(group ~ top_phyloP, scales = 'free_y') + 
  # scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
  xlab('PIP in baselineLF2.2.UKB annotations') + 
  ylab('PIP in Zoonomia annotations') +
  theme_classic(base_size = 6 ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)

pp = ggplot(data = snps_df2, 
            aes(x = PIP_baseline, y = PIP)) + 
  # geom_point(pch = 20, alpha = .25) +
  geom_bin2d(bins = 50) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  # scale_colour_manual(values = phyloP_col) + 
  # guides(colour = guide_legend(nrow = 1)) + 
  facet_grid(group ~ top_phastCons, scales = 'free_y') + 
  # scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
  xlab('PIP in baselineLF2.2.UKB annotations') + 
  ylab('PIP in Zoonomia annotations') +
  theme_classic(base_size = 6 ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)

dev.off()
