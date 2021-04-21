ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)

PROJDIR='code/raw_code/ukbb_gwas'
i_am(file.path(PROJDIR, 'step2_prepare_gwas_pheno.R'))

tab_fn = here('data/tidy_data/tables/readme_ukbb_gwas.tsv')
pheno = read_tsv(tab_fn, col_types = cols()) %>% 
  arrange(FILE) %>% mutate(TRAIT = factor(TRAIT, TRAIT)) %>% arrange(ID)
save_fn = here('data/tidy_data/tables/readme_ukbb_gwas.rds')
saveRDS(pheno, file = save_fn)
