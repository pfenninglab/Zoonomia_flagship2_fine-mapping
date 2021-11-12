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
library(tidymodels); tidymodels_prefer()
library(rmeta)

####################################
## read in table of fine-mapped SNPs ##
PROJDIR='figures/explanatory/finemapping_figures'
i_am(file.path(PROJDIR, 'mainPlots_plot_polyfun_zoonomia_maf.R'))


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
  mutate(inPhyloP = case_when(grepl('Con', top_phyloP) ~ 'PhyloP_Cons.Mam', 
                              grepl('Acc', top_phyloP) ~ 'PhyloP_Acc.Mam',
                              TRUE ~ 'Not in PhyloP'),
         inPhyloP = relevel(factor(inPhyloP), ref = 'Not in PhyloP'),
         inPhastCons = case_when(grepl('Con', top_phastCons) ~ 'PhastCons.Prim', 
                                 TRUE ~ 'Not in PhastCons'),
         inPhastCons = relevel(factor(inPhastCons), ref = 'Not in PhastCons'),
         cCRE_group = relevel(cCRE_group, ref = 'Other'))

table(snps_df$group)

snps_df_mafbin = snps_df %>% 
  select(c(TRAIT:ID, contains('MAF'), zoonomia_CHARs:inPhastCons)) %>%
  pivot_longer(cols = MAFbin_lowfreq_1:MAFbin_frequent_10, 
               names_to = 'MAF_bin', values_to = 'value') %>% 
  filter(value > 0) %>% mutate(
    MAF_bin = factor(MAF_bin, c(paste0('MAFbin_lowfreq_',1:10), paste0('MAFbin_frequent_',1:10)))
  ) %>% select(-value)



##########################################################
## 1. make the main ECDF plots showing PIP improvements ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_fig = 1.75; width_fig = 2.25; font_fig = 7
height_ppt = 4; width_ppt = 8;

## ecdf for phyloP annotations
pp1 = ggplot(snps_df_mafbin %>% filter(group != 'ZoonomiaAnnot', PIP > .1), 
             aes(x = MAF_bin, y = PIP)) + 
  geom_boxplot(aes(fill = group)) + 
  facet_grid(TRAIT~inPhyloP, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_fill_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank())  +
  guides(color = guide_legend(nrow = 2, title.position = 'top'))


## put the eCDF panels together
plot_fn = here(PROJDIR,'plots',
               paste0('polyfun_zoonomia_finemapping_maf_20210727.fig.pdf'))
pdf(plot_fn, height = 12, width = width_ppt, onefile = F)
print(pp1)
dev.off()


