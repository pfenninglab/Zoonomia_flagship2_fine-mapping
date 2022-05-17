ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(data.table)
library(rtracklayer)
library(tidymodels)

####################################
## read in table of fine-mapped SNPs ##
PROJDIR='figures/exploratory/comparing_improved_PIP_scores'
i_am(file.path(PROJDIR, 'step2_check_UNICORNs_overlap.R'))

########################################
## read in the fine-mapping dataframe ##
poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20220517.rds')
snps_df = readRDS(file = poly_fn)

######################
## read in UNICORNs
unicorn_gr = fread(file.path(PROJDIR, 'tables', 'UNICORNs_refined.bed'),
                   header = F, col.names = c('seqnames', 'start', 'end')) %>% GRanges()
names(unicorn_gr) = paste0(seqnames(unicorn_gr), ':', start(unicorn_gr), '-', end(unicorn_gr))

snps_gr = snps_df %>% mutate(seqnames = paste0('chr', CHR), start = POS_hg38, end = POS_hg38) %>%
  select(c(seqnames, start, end)) %>% GRanges()

oo = findOverlaps(subject = snps_gr, query = unicorn_gr)

snps_df2 = snps_df[subjectHits(oo), ] %>%
  mutate(UNICORN = names(unicorn_gr)[queryHits(oo)]) %>% distinct() %>%
  filter(group == 'base + ZooAnnot + cCRE') %>% 
  relocate(UNICORN, .after = SNP) %>% group_by(SNP) %>% 
  mutate(tmp = sum(PIP)) %>% ungroup() %>% 
  arrange(desc(tmp)) %>% select(-tmp) 

snps_df2 %>% writexl::write_xlsx(file.path(PROJDIR, 'tables', 'UNICORN_subset_polyfun_finemapped_snps_zoonomia_20220517.xlsx'))
snps_df2 %>% pull(PIP) %>% mean() #0.1503051
snps_df2 %>% group_by(TRAIT) %>% summarize(score = mean(PIP)) %>% 
  pull(score) %>% summary()


###########################################
## read in random unannotated region subset
nonConUnannotated_gr = fread(file.path(PROJDIR, 'tables', 'unannotated_intergenic_regions_non_constraint_SUBSET.bed'),
                   header = T) %>% rename( 'intergen_chrom' = 'seqnames') %>% GRanges()
names(nonConUnannotated_gr) = paste0(seqnames(nonConUnannotated_gr), ':', start(nonConUnannotated_gr), '-', end(nonConUnannotated_gr))

oo2 = findOverlaps(subject = snps_gr, query = nonConUnannotated_gr)

snps_df3 = snps_df[subjectHits(oo2), ] %>%
  mutate(nonConUnannotated = names(unicorn_gr)[queryHits(oo2)]) %>% distinct() %>%
  filter(group == 'base + ZooAnnot + cCRE') %>% 
  relocate(nonConUnannotated, .after = SNP) %>% group_by(SNP) %>% 
  mutate(tmp = sum(PIP)) %>% ungroup() %>% 
  arrange(desc(tmp)) %>% select(-tmp) 

snps_df3 %>% pull(PIP) %>% mean() #0.05121605
snps_df3 %>% group_by(TRAIT) %>% summarize(score = mean(PIP)) %>% 
  pull(score) %>% summary()

############################################################################
## test PIP of UNICORN snps higher than non-conserved unannotated regions 
snps_df4 = snps_df %>%  
  filter(group == 'base + ZooAnnot + cCRE') %>% 
  mutate(
  UNICORN_group = case_when(SNP %in% snps_df2$SNP ~ "UNICORN",  SNP %in% snps_df3$SNP ~ "Unnannotated",   TRUE ~ 'Other'), 
  UNICORN_group = factor(UNICORN_group, levels = c('UNICORN', "Unnannotated", 'Other')))

snps_df4 %>%
  mutate(name2 = paste(CHR, POS_hg38, SNP, A1, A2, sep = ':')) %>%
  group_by(UNICORN_group) %>%
  filter(!duplicated(name2)) %>%
  summarize(
    numSNP = n(),
    meanPIP = mean(PIP))

snps_df4 %>% nest(data = -c(TRAIT)) %>%
  mutate(
    tmp = map(data,~lm( PIP ~ UNICORN_group, data = .x)),
    tidied = map(tmp, tidy),
    numSNP = map_dbl(data,nrow),
  ) %>% 
  unnest(tidied) %>% arrange(p.value) %>% 
  filter(term != '(Intercept)' ) %>% 
  mutate(FDR = p.adjust(p.value, 'fdr'), 
         term = gsub('UNICORN_group', 'UNICORN vs. ', term)) %>%
  dplyr::select(-c(tmp, data)) %>% 
  as.data.frame() %>% 
  writexl::write_xlsx(file.path(PROJDIR, 'tables', 
                                'UNICORN_vs._nonConstrain_PIP_diffTest_byTrait_20220517.xlsx'))





