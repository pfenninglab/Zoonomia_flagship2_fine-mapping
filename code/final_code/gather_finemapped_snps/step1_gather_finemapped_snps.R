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

CODEDIR='code/final_code/gather_finemapped_snps'
DATADIR='data/tidy_data/polyfun'
PLOTDIR='figures/exploratory/functional_polyfun_Zoonomia_phyloP'
i_am(file.path(CODEDIR, 'step1_gather_finemapped_snps.R'))

###################################
# read in the gwas phenotypes ##
gwas_fn = here('data/tidy_data/tables/readme_ukbb_gwas.rds')
pheno = readRDS(file = gwas_fn)

###############################
# read in annotated PIP SNPs ##
main_groups = c('non-functional','ZoonomiaAnnot','baselineLF2.2.UKB', 'base + ZooAnnot + cCRE')
annot_fn = list.dirs(path = 'data/raw_data', recursive = F) %>%
  lapply(list.dirs, recursive = F) %>% unlist() %>% 
  grep(pattern = 'funct', value = T) %>%
  grep(pattern = 'UKB|annot/|nonfunct|LF2/', value = T) %>%
  lapply(list.files, pattern = '_top_annot.txt.gz', full.names = T) %>%
  unlist()
names(annot_fn) = annot_fn
input = annot_fn %>% lapply(fread) %>% rbindlist(idcol = 'file', fill = TRUE) %>%
  as_tibble() %>% 
  mutate(FILE = file %>% ss('/', 4) %>% ss( '-'),
         value = paste(FILE, SNP, A1, A2, sep = ':'),
         group = ss(file, '/', 3), 
         group = case_when(
           grepl('functional_polyfun_Zoonomia_annot_baselineLF2', group) ~ 'base + ZooAnnot + cCRE', 
           grepl('functional_polyfun_baseline-LF2.2.UKB', group) ~ 'baselineLF2.2.UKB', 
           grepl('functional_polyfun_Zoonomia_annot', group) ~ 'ZoonomiaAnnot', 
           TRUE ~ 'non-functional'),
         group = factor(group, main_groups),
         value = paste(group,FILE, SNP, A1, A2, sep = ':')) %>% select(-file)

# read in causal SNP set ##
annot_fn2 = list.dirs(path = 'data/raw_data', recursive = F) %>%
  lapply(list.dirs, recursive = F) %>% unlist() %>% 
  grep(pattern = 'funct', value = T) %>%
  grep(pattern = 'UKB|annot/|nonfunct|LF2/', value = T) %>%
  lapply(list.files, pattern = 'causal_set.txt.gz', full.names = T) %>%
  unlist()
names(annot_fn2) = annot_fn2
input2 = annot_fn2 %>% lapply(fread) %>% rbindlist(idcol = 'file', fill = TRUE) %>%
  as_tibble() %>% 
  mutate(FILE = file %>% ss('/', 4) %>% ss( '-'),
         group = ss(file, '/', 3), 
         group = case_when(
           grepl('functional_polyfun_Zoonomia_annot_baselineLF2', group) ~ 'base + ZooAnnot + cCRE', 
           grepl('functional_polyfun_baseline-LF2.2.UKB', group) ~ 'baselineLF2.2.UKB', 
           grepl('functional_polyfun_Zoonomia_annot', group) ~ 'ZoonomiaAnnot', 
           TRUE ~ 'non-functional'),
         group = factor(group, main_groups),
         value = paste(group,FILE, SNP, A1, A2, sep = ':')) %>% select(-file)

snps_df = input %>% left_join(input2, by = c('group','FILE', 'SNP','CHR','BP', 'A1', 'A2','PIP','value')) %>% 
  left_join(pheno, by = 'FILE') %>%
  relocate(c(TRAIT,FILE), .before = everything()) %>%
  relocate(Zoonomia_phastCons.43prim.fdr.05:base, .after = everything()) %>%
  relocate(contains('ENCODE3'), .after = contains('Zoonomia')) %>%
  relocate(PIP, .after = P)

snps_df %>% count(group)

########################################
# liftOver hg19 SNP positions to hg38 ##
chain <- file.path("/home/bnphan/resources/liftOver_chainz",
                   'hg19ToHg38.over.chain') %>% import.chain()
snpRanges_hg19 = snps_df %>% 
  mutate(name = paste0('chr',CHR,':',BP)) %>% 
  select(value, name) %>% deframe() %>% GRanges()

snpRanges_hg38 = snpRanges_hg19 %>% rtracklayer::liftOver(chain = chain) %>%
  GenomicRanges::reduce() %>% as.data.frame() %>% rename('group_name' = 'value') %>%
  select(c(value:start))

snps_df = snps_df %>% left_join(snpRanges_hg38, by = 'value') %>%
  rename('start' = 'POS_hg38', 'BP' = 'POS_hg19') %>% 
  filter( seqnames %>% as.character() %>% ss('chr', 2) == CHR) %>% 
  relocate(POS_hg38, .after = 'CHR') %>% select(-c(seqnames, value)) %>%
  mutate(name = paste0(FILE,':',SNP,':chr',CHR,':',POS_hg38,':',A1,':',A2)) %>% 
  relocate(name, .before = CHR) 

#########################################################
## find the finemapped SNPs overlapping CHARs and HARs ##
bed_fn = list.files(path = here('/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/phyloP_cutoff_regions'),
                    pattern = '.bed.gz', full.names = T) %>%
  grep(pattern = 'HAR',value = T) %>% grep(pattern = '20210402|20210416',value = T)
names(bed_fn) = basename(bed_fn) %>% ss('\\_20', 1) 
bed_fn = bed_fn[! grepl('flanking', bed_fn)]
bed_gr = bed_fn %>% lapply(import)

snps_df = bed_gr %>% as.list() %>%
  map( ~ countOverlaps(subject = .x, query = GRanges(paste0('chr',snps_df$CHR,':',snps_df$POS_hg38)))) %>%
  bind_cols(snps_df, . ) %>%
  relocate(all_of(names(bed_gr)), .after = contains('Zoonomia'))

snps_df %>% count(zoonomia_CHARs, group)
snps_df %>% count(zoonomia_HARs, group)

########################################################################
## annotate finemapped SNPs w/ phyloP, primate PhastCons, HARS, CHARs ##
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

snps_df = snps_df %>% mutate(
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
    zoonomia_HARs == 1 ~ 'HAR',
    zoonomia_CHARs == 1 ~ 'CHAR',
    TRUE ~ 'Other'), 
  HAR_group = factor(HAR_group, HAR_lvls)) %>%
  relocate(c(top_phyloP, top_phastCons, cCRE_group, HAR_group), .after = contains('Zoonomia'))

####################################
## save table of fine-mapped SNPs ##
poly_fn = here('data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20220517.rds')
saveRDS(snps_df, file = poly_fn)

snps_df2 = snps_df %>% filter(group =='base + ZooAnnot + cCRE') 
snps_dfList = snps_df2 %>% split(snps_df2$FILE)

snps_tsv = here('data/tidy_data/polyfun/tables/',
                paste0('polyfun_finemapped_snps.',names(snps_dfList),'.20220517.txt.gz'))
tmp = map2(.x = snps_dfList, .y = snps_tsv, ~ write_tsv(.x, file= .y))

#########################################################
## find the finemapped SNPs overlapping CHARs and HARs ##
bed_fn = list.files(path = here('/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/phyloP_cutoff_regions'),
                    pattern = '.bed.gz', full.names = T) %>%
  grep(pattern = 'HAR',value = T) %>% grep(pattern = '20210402|20210416',value = T)
names(bed_fn) = basename(bed_fn) %>% ss('\\_20', 1) 
bed_fn = bed_fn[! grepl('flanking', bed_fn)]
bed_gr = bed_fn %>% lapply(import)

snps_df2 = snps_df %>% dplyr::select(-contains('HAR')) %>% filter(group =='base + ZooAnnot + cCRE')

overLapsHARexact = bed_gr %>%
  map( ~ countOverlaps(subject = .x, query = GRanges(paste0('chr',snps_df2$CHR,':',snps_df2$POS_hg38)))) %>%
  bind_cols(snps_df2, . ) %>%
  mutate(group = 'SNP_in_C/HAR') %>%
  relocate(group, .after = A2)

bed_gr2 = bed_gr %>% lapply(function(gr) {
  start(gr) = round((start(gr) + end(gr)) / 2) 
  end(gr) = start(gr)
  return(gr)
})

overLapsHAR750bp = bed_gr2  %>%
  map( ~ countOverlaps(subject = .x, query = GRanges(paste0('chr',snps_df2$CHR,':',snps_df2$POS_hg38)),
                       maxgap = 750)) %>%
  bind_cols(snps_df2, . ) %>%
  mutate(group = 'SNP750bp_from_C/HAR_center') %>%
  relocate(group, .after = A2)

## combine the two overlapping HARs/CHARs groups
overLapsHAR = overLapsHARexact %>% 
  bind_rows(overLapsHAR750bp) %>% filter_at(vars(contains('HAR')), any_vars(.==1)) %>% 
  mutate(CHAR = paste0('chr',CHR,':',POS_hg38) %>% sapply(function(x){
    tmp = bed_gr2[['zoonomia_CHARs']]
    oo = findOverlaps(subject =  GRanges(x), query = tmp, maxgap = 750)
    return( unique(mcols(tmp)$name[queryHits(oo)])[1])
  }),
  HAR = paste0('chr',CHR,':',POS_hg38) %>% sapply(function(x){
    tmp = bed_gr2[['zoonomia_HARs']]
    oo = findOverlaps(subject =  GRanges(x), query = tmp, maxgap = 750)
    return( unique(mcols(tmp)$name[queryHits(oo)])[1])
  })) %>%
  relocate(c(CHAR, HAR, TRAIT), .after = A2) %>%
  mutate_if(is.character,replace_na,'') %>%
  mutate(group = factor(group, c('SNP_in_C/HAR', 'SNP750bp_from_C/HAR_center'))) %>%
  arrange(group, HAR, CHAR, TRAIT, desc(PIP)) %>%
  distinct(TRAIT, SNP, PIP, .keep_all = T)

with(overLapsHAR %>% mutate(TRAIT = droplevels(TRAIT)), table(TRAIT, group))

#########################################################
## output the HAR/CHAR overlapped SNPs
dir.create(here(DATADIR, 'tables'), showWarnings = F)
har_tsv = here(DATADIR, 'tables',
               'polyfun_UKBB_finemapped_snps_overlap_HAR_CHAR_20220517.tsv')
write_tsv(overLapsHAR, file = har_tsv)

dir.create(here(DATADIR, 'rdas'), showWarnings = F)
har_rds = here(DATADIR, 'rdas',
               'polyfun_UKBB_finemapped_snps_overlap_HAR_CHAR_20220517.rds')
saveRDS(overLapsHAR, file = har_rds)
