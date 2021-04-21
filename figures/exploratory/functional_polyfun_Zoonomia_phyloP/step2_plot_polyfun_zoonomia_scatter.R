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

poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia.rds')
snps_df = readRDS(file = poly_fn)

snps_df2 = snps_df %>%
  group_by(TRAIT, SNP) %>% filter(n() == 4) %>%
  mutate(PIP_baseline = PIP[which(group == 'baselineLF2.2.UKB')]) %>% 
  ungroup() %>% filter(PIP_baseline > 0.01) %>% 
  filter(! group %in% c( 'baselineLF2.2.UKB'))

#################################
## make plots for presentation ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
phyloP_col = c(brewer.pal(4,'Blues'),brewer.pal(4,'Reds'),'#33a02c')
names(phyloP_col) = c('Cons_0.10','Cons_0.05', 'Cons_0.01', 'Cons_0.001', 
                     'Accl_0.10','Accl_0.05','Accl_0.01','Accl_0.001','Other')


plot_fn = here(PROJDIR,'plots',
               paste0('polyfun_zoonomia_finemapping_PIPscatter_20210405.ppt.pdf'))
pdf(plot_fn, height = height_ppt, width = width_ppt)
pp = ggplot(data = snps_df2, 
            aes(x = PIP_baseline, y = PIP)) + 
  # geom_point(pch = 20, alpha = .25) +
  geom_bin2d(bins = 50) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  # scale_colour_manual(values = phyloP_col) + 
  # guides(colour = guide_legend(nrow = 1)) + 
  facet_grid(group ~ phyloP, scales = 'free_y') + 
  # scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
  xlab('PIP in baselineLF2.2.UKB annotations') + 
  ylab('PIP in Zoonomia annotations') +
  theme_classic(base_size = 6 ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)

pp = ggplot(data = snps_df2 %>% 
              mutate(phyloP = ifelse(phyloP == 'Accl_0.10', 'Other', as.character(phyloP)),
                     phyloP = factor(phyloP, levels(snps_df$phyloP))), 
            aes(x = PIP_baseline, y = PIP)) + 
  # geom_point(pch = 20, alpha = .25) +
  geom_bin2d(bins = 50) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  # scale_colour_manual(values = phyloP_col) + 
  # guides(colour = guide_legend(nrow = 1)) + 
  facet_grid(group ~ phyloP, scales = 'free_y') + 
  # scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') +
  xlab('PIP in baselineLF2.2.UKB annotations') + 
  ylab('PIP in Zoonomia annotations') +
  theme_classic(base_size = 6 ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp)
dev.off()


#################################
## make plots for presentation ##
snps_df2 = snps_df %>% count(FILE, phyloP)
snps_df %>% filter(inHAR != 'notInHAR') %>% pull(SNP) 
