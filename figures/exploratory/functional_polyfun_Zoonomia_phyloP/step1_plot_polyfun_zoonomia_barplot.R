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


####################################
# read in annotated 0.95 PIP SNPs ##
if(FALSE){
  gwas_fn = here('data/tidy_data/tables/readme_ukbb_gwas.rds')
  pheno = readRDS(file = gwas_fn)
  
  annot_fn = list.dirs(path = 'data/raw_data', recursive = F) %>%
    lapply(list.dirs, recursive = F) %>% unlist() %>% 
    grep(pattern = 'funct', value = T) %>%
    grep(pattern = 'UKB|annot/|nonfunct|LF2/', value = T) %>%
    lapply(list.files, pattern = '_top_annot.txt.gz', full.names = T) %>%
    unlist()
  names(annot_fn) = annot_fn
  input = annot_fn %>% lapply(fread) %>% rbindlist(idcol = 'file', fill = TRUE) %>%
    as_tibble() 
  
  main_groups = c('non-functional','ZoonomiaAnnot','baselineLF2.2.UKB', 'base + ZooAnnot + cCRE')
  snps_df = input %>% 
    select(c(file:A2, contains('Zoonomia'),contains("ENCODE3") ,contains('synonymous'))) %>%
    mutate(FILE = ss(ss(file, '/', 4), '-'),
           group = ss(file, '/', 3), 
           group = case_when(
             grepl('functional_polyfun_Zoonomia_annot_baselineLF2', group) ~ 'base + ZooAnnot + cCRE', 
             grepl('functional_polyfun_baseline-LF2.2.UKB', group) ~ 'baselineLF2.2.UKB', 
             grepl('functional_polyfun_Zoonomia_annot', group) ~ 'ZoonomiaAnnot', 
             TRUE ~ 'non-functional'),
           group = factor(group, main_groups)) %>%
    # filter(group != 'base + ZooAnnot + cCRE') %>%
    left_join(pheno, by = 'FILE')
  
  table(snps_df$group)
  
  ######################################################################
  ## read in the bed files for phyloP, primate PhastCons, HARS, CHARs ##
  bed_fn = list.files(path = here('data/raw_data/zoonomia_annotations/bed_hg19'),
                      pattern = '.bed.gz', full.names = T) %>%
    grep(pattern = 'HAR',value = T)
  names(bed_fn) = basename(bed_fn) %>% ss('\\.hg19\\.bed\\.gz') 
  bed_fn = bed_fn[! grepl('flanking', bed_fn)]
  bed_gr = bed_fn %>% lapply(import) %>% GRangesList()
  
  snps_df = bed_gr %>% as.list() %>%
    map( ~ countOverlaps(subject = .x, query = GRanges(paste0('chr',snps_df$CHR,':',snps_df$BP)))) %>%
    bind_cols(snps_df, . )
  
  ########################################################################
  ## annotate finemapped SNPs w/ phyloP, primate PhastCons, HARS, CHARs ##
  top_phyloP_lvls = c('Con.top.0-1%', 'Con.top.1-2%', 'Con.top.2-3%', 'Con.top.3-4%',
                      'Acc.top.3-4%', 'Acc.top.2-3%', 'Acc.top.1-2%', 'Acc.top.0-1%',  
                      'Other')
  HAR_lvls = c('HAR1500bp', 'CHAR1500bp', 'Other')
  HAR_cols = c(brewer.pal(3,'Accent'))
  names(HAR_cols) = HAR_lvls
  
  snps_df2 = snps_df %>% mutate(
    # Group 241mammals phyloP
    top_phyloP = case_when(
      `Zoonomia_phyloPaccl.241mam.top0-1%` == 1 ~ 'Acc.top.0-1%',
      `Zoonomia_phyloPaccl.241mam.top1-2%` == 1 ~ 'Acc.top.1-2%',
      `Zoonomia_phyloPaccl.241mam.top2-3%` == 1 ~ 'Acc.top.2-3%',
      `Zoonomia_phyloPaccl.241mam.top3-4%` == 1 ~ 'Acc.top.3-4%',
      `Zoonomia_phyloPcons.241mam.top0-1%` == 1 ~ 'Con.top.0-1%',
      `Zoonomia_phyloPcons.241mam.top1-2%` == 1 ~ 'Con.top.1-2%',
      `Zoonomia_phyloPcons.241mam.top2-3%` == 1 ~ 'Con.top.2-3%',
      `Zoonomia_phyloPcons.241mam.top3-4%` == 1 ~ 'Con.top.3-4%',
      TRUE ~ 'Other'),
    top_phyloP = factor(top_phyloP, top_phyloP_lvls), 
    # Group primate PhastCons
    top_phastCons = case_when(
      `Zoonomia_phastCons.43prim.top0-1%` == 1 ~ 'PhastCons.top.0-1%',
      `Zoonomia_phastCons.43prim.top1-2%` == 1 ~ 'PhastCons.top.1-2%',
      `Zoonomia_phastCons.43prim.top2-3%` == 1 ~ 'PhastCons.top.2-3%',
      `Zoonomia_phastCons.43prim.top3-4%` == 1 ~ 'PhastCons.top.3-4%',
      `Zoonomia_phastCons.43prim.top4-5%` == 1 ~ 'PhastCons.top.4-5%',
      TRUE ~ 'Other'),
    top_phastCons = factor(top_phastCons, top_phastCons_lvls), 
    # Group ENCODE3 cCREs
    cCRE_group = case_when(
      ENCODE3.dELS == 1 ~ 'dELS',
      ENCODE3.pELS == 1 ~ 'pELS',
      ENCODE3.PLS == 1 ~ 'PLS',
      TRUE ~ 'Other'), 
    cCRE_group = factor(cCRE_group, ENCODE3_cCRE_lvls), 
    # Group HAR1500bp or CHAR1500bp not in HAR1500bp
    HAR_group = case_when(
      zoonomia_HARs_20210402.1500bp == 1 ~ 'HAR1500bp',
      zoonomia_CHARs_20210416.1500bp == 1 ~ 'CHAR1500bp',
      TRUE ~ 'Other'), 
    HAR_group = factor(HAR_group, HAR_lvls))
  
  ####################################
  ## save table of fine-mapped SNPs ##
  poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20210520.rds')
  saveRDS(snps_df2, file = poly_fn)
  snps_df = snps_df2
} else{
  poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20210520.rds')
  snps_df = readRDS(file = poly_fn)
}



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
