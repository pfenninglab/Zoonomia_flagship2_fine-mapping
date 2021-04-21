ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)

PROJDIR='figures/exploratory/functional_polyfun_Zoonomia_phyloP'
i_am(file.path(PROJDIR, 'step1_plot_polyfun_zoonomia.R'))

####################################
# read in annotated 0.95 PIP SNPs ##
annot_fn = list.dirs(path = 'data/raw_data', recursive = F) %>%
  lapply(list.dirs, recursive = F) %>% unlist() %>% 
  grep(pattern = 'funct', value = T) %>%
  lapply(list.files, pattern = 'causal_set.txt.gz', full.names = T) %>%
  unlist()
names(annot_fn) = annot_fn

input = annot_fn %>% lapply(fread) %>% rbindlist(idcol = 'file', fill = TRUE) %>%
  as_tibble() 

snps_df = input %>% mutate( 
  group = ss(file, '/', 3), 
  group = case_when(
    grepl('functional_polyfun_Zoonomia_phyloP_phast_HAR', group) ~ 'ZooPhyloPhastHAR', 
    grepl('functional_polyfun_Zoonomia_phyloP', group) ~ 'ZooPhyloHAR', 
    grepl('functional_polyfun_merged_baseline_zoonomia',group) ~ 'base + ZooExt500',
    grepl('functional_polyfun_baseline-LF2.2.UKB', group) ~ 'baselineLF', 
    grepl('functional_polyfun_Zoonomia_conservation', group) ~ 'ZooPPHExt500', 
    grepl('phyloP', group) ~ 'Zoonomia', TRUE ~ 'non-functional'),
  group = factor(group, c('non-functional','ZooPhyloHAR', 'ZooPhyloPhastHAR','ZooPPHExt500',
                          'base + ZooExt500','baselineLF')),
  gwas = ss(ss(file, '/', 4), '-')) %>% group_by(gwas) %>% 
  filter(length(unique(group)) >= 3) %>% ungroup()

#################################
## make plots for presentation ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
annot_col = c(brewer.pal(4,'Blues'),brewer.pal(4,'Reds'),'#33a02c')
names(annot_col) = c('Cons_0.10','Cons_0.05', 'Cons_0.01', 'Cons_0.001', 
                      'Accl_0.10', 'Accl_0.05', 'Accl_0.01', 'Accl_0.001',
                    'Other')

plot_fn = here(PROJDIR,'plots','polyfun_zoonomia_finemapping_20210401.ppt.pdf')
pdf(plot_fn, height = height_ppt, width = width_ppt)
ggplot(data = snps_df , aes(x = group, fill = group)) +
  # geom_bar(aes(fill = annot)) + 
  geom_bar() + 
  geom_text(stat='count', aes(label=..count..), vjust=-1, size= 2.5)+
  # scale_fill_manual(values = annot_col) + 
  scale_fill_carto_d(palette = "Safe") + 
  facet_wrap(~gwas, scales = 'free_y', ncol = 6) + 
  scale_y_continuous(expand = expansion(mult = c(0, .4))) + 
  xlab('Fine-mapping group') + ylab('Number of SNPs w/ PIP >0.95') +
  theme_classic(base_size = 6.5 ) + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), 
        legend.position = 'bottom') 
dev.off()

