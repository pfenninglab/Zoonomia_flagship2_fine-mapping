ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(data.table)


####################################
## read in table of fine-mapped SNPs ##
PROJDIR='figures/exploratory/comparing_improved_PIP_scores'
i_am(file.path(PROJDIR, 'step1_check_snps.R'))

########################################
## read in the fine-mapping dataframe ##
poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20210520.rds')
snps_df = readRDS(file = poly_fn)



snps_df %>% filter(SNP == 'rs76488452', TRAIT== 'BMI') %>% pull(P)


## Sent to Steven Gazal 8/31/21
## I grabbed the SNPs from LDproxy with that SNP from here. It looks like that SNP is weakly in LD with many of the P< 10^-14 SNPs (R2 ~ 0.10- 0.20). SuSIE and any fine-mapper would have accounted for the weak LD and probalby the annotatiosn for this SNP is really strong for the heritability. Does that help?
sumstats_bmi = fread(here('data/tidy_data/sumstats/body_BMIz.sumstats.gz'))
sumstats_bmi %>% filter(SNP == 'rs76488452')

SNPs_in_LD = fread(file.path(PROJDIR, 'SNP_inLD_with_rs76488452.txt')) %>%
  filter(R2 > .1) %>% rename('SNP' = 'RS_Number','R2_with_rs76488452' = 'R2')

sumstats_bmi %>% inner_join(SNPs_in_LD) %>% arrange(P) %>% 
  write_tsv(file.path(PROJDIR, 'BMI_finemap_snp_rs76488452_and_LD_friends.tsv'))

fread(file.path(PROJDIR, 'SNP_inLD_with_rs76488452.txt')) %>%
  filter(RS_Number %in% c('rs58907687', 'rs76764809'))

############################################################################################################
# Notably, PolyFun with the baseline-LF+Zoonomia model detected 2542 variants constrained in mammals fine-mapped with high confidence (PIP > 0.95) across all the UK Biobank traits, against 2049 and 2486 when using the non-functional and baseline-LF and models, respectively (Figure BNP1 - SUPP). 

snps_df %>% filter(PIP > 0.95, group != 'ZoonomiaAnnot') %>%
  group_by(group) %>% 
  summarise(ConsInX = sum(Zoonomia_phyloPcons.241mam.fdr.05==1), X = n(),
            PercentConsFineMapped = ConsInX/X * 100 )
# group                  ConsInX     X PercentConsFineMapped
# 1 non-functional             598  2048                  29.2
# 2 baselineLF2.2.UKB          833  2486                  33.5
# 3 base + ZooAnnot + cCRE     913  2542                  35.9


snps_df %>% filter(PIP > 0.75, group != 'ZoonomiaAnnot') %>%
  group_by(group) %>% 
  summarise(ConsInX = sum(Zoonomia_phyloPcons.241mam.fdr.05==1), X = n(),
            PercentConsFineMapped = ConsInX/X * 100 )
# group                  ConsInX     X PercentConsFineMapped
# 1 non-functional             732  2988                  24.5
# 2 baselineLF2.2.UKB         1216  4099                  29.7
# 3 base + ZooAnnot + cCRE    1407  4289                  32.8


table_out_fn = here(PROJDIR, 'table_fine-mapping_perTrait_3Models_2021.11.11.xlsx')
snps_df %>% filter(PIP > 0.75, group != 'ZoonomiaAnnot') %>%
  rename('N' = 'N.y') %>%
  group_by(group, TRAIT, N, Neff) %>% 
  summarise(numSNP_PIP_gt_0.75 = n(),
            numConservedSNP = sum(Zoonomia_phyloPcons.241mam.fdr.05==1), 
            PercentConsFineMapped = numConservedSNP / numSNP_PIP_gt_0.75 * 100 ,
            PercentConsFineMapped = signif(PercentConsFineMapped, 2)) %>%
  writexl::write_xlsx(table_out_fn)


snps_df %>% filter(PIP > 0.5, group != 'ZoonomiaAnnot') %>%
  group_by(group) %>% 
  summarise(ConsInX = sum(Zoonomia_phyloPcons.241mam.fdr.05==1), X = n(),
            PercentConsFineMapped = ConsInX/X * 100 )

# group                  ConsInX     X PercentConsFineMapped
# 1 non-functional             916  4564                  20.1
# 2 baselineLF2.2.UKB         1681  6624                  25.4
# 3 base + ZooAnnot + cCRE    1986  6891                  28.8




