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
poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20220517.rds')
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
# <fct>                    <int> <int>                 <dbl>
# 1 non-functional             655  2238                  29.3
# 2 baselineLF2.2.UKB          917  2717                  33.8
# 3 base + ZooAnnot + cCRE    1029  2824                  36.4

snps_df %>% filter(PIP > 0.75, group != 'ZoonomiaAnnot') %>%
  group_by(group) %>% 
  summarise(ConsInX = sum(Zoonomia_phyloPcons.241mam.fdr.05==1), X = n(),
            PercentConsFineMapped = ConsInX/X * 100 )
# group                  ConsInX     X PercentConsFineMapped
# 1 non-functional             823  3318                  24.8
# 2 baselineLF2.2.UKB         1350  4497                  30.0
# 3 base + ZooAnnot + cCRE    1586  4792                  33.1


table_out_fn = here(PROJDIR, 'table_fine-mapping_perTrait_3Models_20220517.xlsx')
snps_df %>% filter(PIP > 0.75, group != 'ZoonomiaAnnot') %>%
  dplyr::rename('N' = 'N.y') %>%
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
# 1 non-functional            1031  5068                  20.3
# 2 baselineLF2.2.UKB         1869  7338                  25.5
# 3 base + ZooAnnot + cCRE    2256  7745                  29.1




