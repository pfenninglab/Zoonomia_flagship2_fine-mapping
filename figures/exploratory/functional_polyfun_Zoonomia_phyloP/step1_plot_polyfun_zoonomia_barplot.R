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

PROJDIR='figures/exploratory/functional_polyfun_Zoonomia_phyloP'
i_am(file.path(PROJDIR, 'step1_plot_polyfun_zoonomia.R'))

###################################
# read in the gwas phenotypes ##
gwas_fn = here('data/tidy_data/tables/readme_ukbb_gwas.rds')
pheno = readRDS(file = gwas_fn)

####################################
# read in annotated 0.95 PIP SNPs ##
annot_fn = list.dirs(path = 'data/raw_data', recursive = F) %>%
  lapply(list.dirs, recursive = F) %>% unlist() %>% 
  grep(pattern = 'funct', value = T) %>%
  lapply(list.files, pattern = '_top_annot.txt.gz', full.names = T) %>%
  unlist()
names(annot_fn) = annot_fn
input = annot_fn %>% lapply(fread) %>% rbindlist(idcol = 'file', fill = TRUE) %>%
  as_tibble() 

main_groups = c('non-functional','ZoonomiaExt500','baselineLF2.2.UKB', 'base + ZooExt500')
snps_df = input %>% select(file:Zoonomia_HAR.extend500) %>% 
  mutate(FILE = ss(ss(file, '/', 4), '-'),
  group = ss(file, '/', 3), 
  group = case_when(
    grepl('functional_polyfun_Zoonomia_phyloP_phast_HAR', group) ~ 'ZooPhylopPhastHAR', 
    grepl('functional_polyfun_Zoonomia_phyloP', group) ~ 'ZooPhylopHAR', 
    grepl('functional_polyfun_merged_baseline_zoonomia',group) ~ 'base + ZooExt500',
    grepl('functional_polyfun_baseline-LF2.2.UKB', group) ~ 'baselineLF2.2.UKB', 
    grepl('functional_polyfun_Zoonomia_conservation', group) ~ 'ZoonomiaExt500', 
    grepl('phyloP', group) ~ 'Zoonomia', TRUE ~ 'non-functional'),
  group = factor(group, c('non-functional','ZooPhylopHAR', 'ZooPhylopPhastHAR', 
                          'ZoonomiaExt500','baselineLF2.2.UKB','base + ZooExt500')),
  # annotations by PhyloP
  phyloP = case_when(
    `Zoonomia_phyloPam_cons.fdr.001` ==1 ~ 'Cons_0.001',
    `Zoonomia_phyloPam_cons.fdr.01`  ==1 ~ 'Cons_0.01',
    `Zoonomia_phyloPam_cons.fdr.05`  ==1 ~ 'Cons_0.05',
    `Zoonomia_phyloPam_cons.fdr.10`  ==1 ~ 'Cons_0.10',
    `Zoonomia_phyloPam_accl.fdr.001` ==1 ~ 'Accl_0.001',
    `Zoonomia_phyloPam_accl.fdr.01`  ==1 ~ 'Accl_0.01',
    `Zoonomia_phyloPam_accl.fdr.05`  ==1 ~ 'Accl_0.05',
    `Zoonomia_phyloPam_accl.fdr.10`  ==1 ~ 'Accl_0.10', TRUE ~ 'Other'
  ) %>% factor(c('Cons_0.001','Cons_0.01', 'Cons_0.05','Cons_0.10',
                 'Accl_0.10', 'Accl_0.05', 'Accl_0.01', 'Accl_0.001','Other')),
  # annotations by PhastCons
  phastCons = case_when(
    `Zoonomia_PhastCons_43primates` == 1 && `Zoonomia_PhastCons_241mammals` == 0 ~ '43prim only',
    `Zoonomia_PhastCons_43primates` == 0 && `Zoonomia_PhastCons_241mammals` == 1 ~ '241mam only',
    `Zoonomia_PhastCons_43primates` == 1 && `Zoonomia_PhastCons_241mammals` == 1 ~ '241mam + 43prim',
    TRUE ~ 'Other') %>% factor(c('43prim only', '241mam only', '241mam + 43prim', 'Other')),
  # annotations by HAR
  inHAR = case_when( `Zoonomia_HAR` ==1 ~ 'HAR',
    `Zoonomia_HAR.extend500` == 1 ~ 'HAR_ext500', TRUE ~ 'notInHAR') %>% 
  factor(c('HAR', 'HAR_ext500', 'notInHAR'))) %>% group_by(FILE) %>% 
  filter(length(unique(group)) >= 3) %>% ungroup() %>% 
  left_join(pheno, by = 'FILE') %>% filter(group %in% main_groups) 

table(snps_df$group, snps_df$phyloP)


####################################
## save table of fine-mapped SNPs ##
poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia.rds')
saveRDS(snps_df, file = poly_fn)


#################################
## make plots for presentation ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
phyloP_col = c(brewer.pal(4,'Blues'),brewer.pal(4,'Reds'),'#33a02c')
names(phyloP_col) = c('Cons_0.10','Cons_0.05', 'Cons_0.01', 'Cons_0.001', 
                     'Accl_0.10','Accl_0.05','Accl_0.01','Accl_0.001','Other')


plot_fn = here(PROJDIR,'plots',
               paste0('polyfun_zoonomia_finemapping_20210402.ppt.pdf'))
pdf(plot_fn, height = height_ppt, width = width_ppt)
for (cutoff in c(.95, .9, .5, .25)){
plot_fn = here(PROJDIR,'plots',
  paste0('polyfun_zoonomia_finemapping_PIP',cutoff,'_20210402.ppt.pdf'))
pdf(plot_fn, height = height_ppt, width = width_ppt)
pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = group)) +
  geom_bar(aes(fill = phyloP)) + 
  geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
  scale_fill_manual(values = phyloP_col) + 
  facet_wrap(~TRAIT, scales = 'free_y', ncol = 6) + 
  scale_y_continuous(expand = expansion(mult = c(0, .4))) + 
  xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
  theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
print(pp); dev.off(); print(pp)
}
dev.off()


#################################
## make plots for presentation ##
snps_df2 = snps_df %>% count(FILE, phyloP)
snps_df %>% filter(inHAR != 'notInHAR') %>% pull(SNP) 
