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
i_am(file.path(PROJDIR, 'step3_plot_polyfun_zoonomia_ecdf.R'))


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


##########################################################
## 1. make the main ECDF plots showing PIP improvements ##
dir.create(here(file.path(PROJDIR, 'plots')), showWarnings = F)
height_fig = 1.75; width_fig = 2.25; font_fig = 7
height_ppt = 4; width_ppt = 8;

## ecdf for phyloP annotations
pp1 = ggplot(snps_df %>% filter(group != 'ZoonomiaAnnot'),
             aes(x = PIP, color = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(~inPhyloP, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_colour_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank())  +
  guides(color = guide_legend(nrow = 2, title.position = 'top'))

## ecdf for phastCons annotations
pp2 = ggplot(snps_df %>% filter(group != 'ZoonomiaAnnot'),
            aes(x = PIP, color = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(~inPhastCons, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_colour_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, title.position = 'top'))

## ecdf for cCRE annotations
pp3 = ggplot(snps_df %>% filter(group != 'ZoonomiaAnnot'),
             aes(x = PIP, color = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(~cCRE_group, scales = 'free_y') + 
  ylab('CDF Proportion') + ylim(c(.8, NA)) + 
  scale_colour_carto_d(palette = "Vivid") +
  theme_classic(base_size = font_fig ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.key.size = unit(.25, 'cm'), legend.position = 'bottom',
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,0),
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, title.position = 'top'))

## put the eCDF panels together
plot_fn = here(PROJDIR,'plots',
               paste0('polyfun_zoonomia_finemapping_ecdf_20210724.fig.pdf'))
pdf(plot_fn, height = 4, width = width_fig, onefile = F)
pp = ggarrange(pp1, pp2, pp3, 
          font.label = list(size = font_fig + 4, color = "black", face = "bold"),
          common.legend = TRUE, legend = 'top', ncol = 1, nrow = 3)
print(pp)
dev.off()


####################
## ecdf KS tests ##
ksTests = snps_df %>% filter(grepl('base', group)) %>%
  nest(data = -c(TRAIT)) %>%
  mutate(fit = map(data, ~ ks.test(x = .x %>% filter(group == 'baselineLF2.2.UKB') %>% pull(PIP), 
                                   y = .x %>% filter(group == 'base + ZooAnnot + cCRE') %>% pull(PIP), 
                                   alternative = 'greater')),
         tidied = map(fit, tidy)) %>% 
  unnest(tidied) %>%  select(-data, -fit) %>%  
  # mutate(group = 'All_SNPs') %>% 
  relocate(group, .after = TRAIT) %>%
  arrange(p.value)

ksTests2 = snps_df %>% filter(grepl('base', group)) %>%
  nest(data = -c(inPhyloP)) %>%
  mutate(fit = map(data, ~ ks.test(x = .x %>% filter(group == 'baselineLF2.2.UKB') %>% pull(PIP), 
                                   y = .x %>% filter(group == 'base + ZooAnnot + cCRE') %>% pull(PIP), 
                                   alternative = 'greater')),
         tidied = map(fit, tidy)) %>% 
  unnest(tidied) %>%  select(-data, -fit) %>%  
  rename(`group` = `inPhyloP`) %>%
  # mutate(FDR = p.adjust(p.value, 'fdr')) %>%
  arrange(p.value)

ksTests3 = snps_df %>% filter(grepl('base', group)) %>%
  nest(data = -c(inPhastCons)) %>%
  mutate(fit = map(data, ~ ks.test(x = .x %>% filter(group == 'baselineLF2.2.UKB') %>% pull(PIP), 
                                   y = .x %>% filter(group == 'base + ZooAnnot + cCRE') %>% pull(PIP), 
                                   alternative = 'greater')),
         tidied = map(fit, tidy)) %>% 
  unnest(tidied) %>%  select(-data, -fit) %>%  
  rename(`group` = `inPhastCons`) %>%
  # mutate(FDR = p.adjust(p.value, 'fdr')) %>%
  arrange(p.value)

ksTests4 = snps_df %>% filter(grepl('base', group)) %>%
  nest(data = -c(cCRE_group)) %>%
  mutate(fit = map(data, ~ ks.test(x = .x %>% filter(group == 'baselineLF2.2.UKB') %>% pull(PIP), 
                                   y = .x %>% filter(group == 'base + ZooAnnot + cCRE') %>% pull(PIP), 
                                   alternative = 'greater')),
         tidied = map(fit, tidy)) %>% 
  unnest(tidied) %>%  select(-data, -fit) %>% 
  rename(`group` = `cCRE_group`) %>%
  # mutate(FDR = p.adjust(p.value, 'fdr')) %>%
  arrange(p.value)

ksTests %>% writexl::write_xlsx( here(PROJDIR,'tables', 'ksTest.PIP.byTrait.results.xlsx'))
ksTests2 %>% writexl::write_xlsx( here(PROJDIR,'tables', 'ksTest.PIP.byTraitAndPhyloP.results.xlsx'))
ksTests3 %>% writexl::write_xlsx( here(PROJDIR,'tables', 'ksTest.PIP.byTraitAndPhastCons.results.xlsx'))

rbind( ksTests2, ksTests3, ksTests4)%>% 
  mutate(p.bonferroni = p.adjust(p.value, 'bonferroni'),
         signif.sign = case_when(
           p.bonferroni < 0.0001 ~ '***',
           p.bonferroni < 0.001 ~ '**',
           p.bonferroni < 0.01 ~ '*',
           TRUE~ 'NS' )) %>%
  writexl::write_xlsx( here(PROJDIR,'tables', 'Table_BNP2-1_ksTest.results.PIP.byAnnotations.xlsx'))

###########################################
## 3. count table of SNPs in categories ###
dat1 = snps_df %>% filter(PIP > .1) %>%
  count(TRAIT, inPhyloP, group,  name = 'num') %>%
  group_by(group) %>% mutate(freq = 100 * num/sum(num))

dat2 = snps_df %>% filter(PIP > .1) %>%
  count(TRAIT, inPhastCons, group,  name = 'num') %>%
  group_by(group) %>% mutate(freq = 100 * num/sum(num))

dat3 = snps_df %>% filter(PIP > .1) %>%
  count(TRAIT, cCRE_group, group,  name = 'num') %>%
  group_by(group) %>% mutate(freq = 100 * num/sum(num))

dat4 = snps_df %>%
  count(TRAIT, inPhyloP, group,  name = 'num') %>%
  group_by(group) %>% mutate(freq = 100 * num/sum(num))

dat5 = snps_df %>%
  count(TRAIT, inPhastCons, group,  name = 'num') %>%
  group_by(group) %>% mutate(freq = 100 * num/sum(num))

dat6 = snps_df %>%
  count(TRAIT, cCRE_group, group,  name = 'num') %>%
  group_by(group) %>% mutate(freq = 100 * num/sum(num))

datList = list(dat1, )

###################################
## plot of number of SNPs more ###
label_cols = c('PhyloP_Cons.Mam' = "darkblue", 'PhyloP_Acc.Mam' = 'darkred', 
               'PhastCons.Prim' = 'lightgreen', 'PLS' = 'bisque', 
               'pELS' = 'purple', 'dELS' = 'coral')

pp1 = ggplot(dat1 %>% select(-freq) %>% spread(group, num) %>% 
               filter(inPhyloP != 'Not in PhyloP'), 
       aes(x =  `baselineLF2.2.UKB`, y =  `base + ZooAnnot + cCRE`)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') + 
  geom_point(pch = 16, aes(color = inPhyloP)) + 
  theme_classic(base_size = font_fig ) + 
  ylab('base + ZooAnnot  ') + xlab('baselineLF v2.2 UKB') + 
  ggtitle('# SNPs, PIP > 0.1') + 
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = label_cols) +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(.04, 'cm'),         
        legend.key.width=unit(.5,"line")) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))


pp2 = ggplot(dat2 %>% select(-freq) %>% spread(group, num) %>% 
               filter(inPhastCons != 'Not in PhastCons'), 
             aes(x =  `baselineLF2.2.UKB`, y =  `base + ZooAnnot + cCRE`)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') + 
  geom_point(pch = 16, aes(colour = 'PhastCons.Prim')) + 
  theme_classic(base_size = font_fig ) + 
  ylab('base + ZooAnnot        ') + xlab('baselineLF v2.2 UKB') + 
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = label_cols) +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), legend.title = element_blank(), 
        legend.spacing.x = unit(.04, 'cm'), legend.key.width=unit(.5,"line")) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))


pp3 = ggplot(dat3 %>% select(-freq) %>% spread(group, num) %>% 
               filter(cCRE_group != 'Other'), 
             aes(x =  `baselineLF2.2.UKB`, y =  `base + ZooAnnot + cCRE`)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') + 
  geom_point(pch = 16, aes(color = cCRE_group)) + 
  theme_classic(base_size = font_fig ) + 
  ylab('base + ZooAnnot        ') + xlab('baselineLF v2.2 UKB') + 
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = label_cols) +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), legend.title = element_blank(), 
        legend.spacing.x = unit(.04, 'cm'), legend.key.width=unit(.5,"line")) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))


##################################################
## plot of perentage of SNPs finemapped better ###
pp4 = ggplot(dat4 %>% select(-num) %>% spread(group, freq) %>% 
               filter(inPhyloP != 'Not in PhyloP'), 
             aes(x =  100 * (`baselineLF2.2.UKB` - `non-functional`) / `non-functional`,
                 y =  100 * (`base + ZooAnnot + cCRE` -  `non-functional`)  / `non-functional`)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') + 
  geom_point(pch = 16, aes(color = inPhyloP)) + 
  theme_classic(base_size = font_fig ) + 
  ylab('base + ZooAnnot  ') + xlab('baselineLF v2.2 UKB') + 
  ggtitle('% more vs. nonfunctional') + 
  scale_color_manual(values = label_cols) +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.spacing.x = unit(.04, 'cm'), 
        legend.key.width=unit(.5,"line")) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))


pp5 = ggplot(dat5 %>% select(-num) %>% spread(group, freq) %>% 
               filter(inPhastCons != 'Not in PhastCons'), 
             aes(x =  100 * (`baselineLF2.2.UKB` - `non-functional`) / `non-functional`,
                 y =  100 * (`base + ZooAnnot + cCRE` -  `non-functional`)  / `non-functional`)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') + 
  geom_point(pch = 16 , aes(colour = 'PhastCons.Prim')) + 
  theme_classic(base_size = font_fig ) + 
  ylab('base + ZooAnnot        ') + xlab('baselineLF v2.2 UKB') + 
  scale_color_manual(values = label_cols) +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), legend.title = element_blank(), 
        legend.spacing.x = unit(.04, 'cm'), legend.key.width=unit(.5,"line")) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))


pp6 = ggplot(dat6 %>% select(-num) %>% spread(group, freq) %>% 
               filter(cCRE_group != 'Other'), 
             aes(x =  100 * (`baselineLF2.2.UKB` - `non-functional`) / `non-functional`,
                 y =  100 * (`base + ZooAnnot + cCRE` -  `non-functional`)  / `non-functional`)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'red') + 
  geom_point(pch = 16, aes(color = cCRE_group)) + 
  theme_classic(base_size = font_fig ) + 
  ylab('base + ZooAnnot        ') + xlab('baselineLF v2.2 UKB') + 
  scale_color_manual(values = label_cols) +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5), legend.title = element_blank(), 
        legend.spacing.x = unit(.04, 'cm'), legend.key.width=unit(.5,"line")) +
  guides(color = guide_legend(nrow = 1, title.position = 'top'))

# put together the plot panels, gain in fine-mapped SNPs
plot_fn3 = here(PROJDIR,'plots',
                paste0('polyfun_zoonomia_finemapping_baseVsCombine_20210724.fig.pdf'))
pdf(plot_fn3, height = 4, width = 3.5, onefile = F)

pa = ggarrange(pp1,pp4, legend = 'bottom', ncol = 2, nrow = 1, 
               font.label = list(size = font_fig + 4, color = "black", face = "bold"),
               common.legend = TRUE)


pb = ggarrange(pp2,pp5,legend = 'bottom', ncol = 2, nrow = 1,
               font.label = list(size = font_fig + 4, color = "black", face = "bold"),
               common.legend = TRUE)


pc = ggarrange(pp3, pp6,legend = 'bottom', ncol = 2, nrow = 1,
               font.label = list(size = font_fig + 4, color = "black", face = "bold"),
               common.legend = TRUE)
ggarrange(pa, pb, pc, ncol = 1, nrow = 3,legend = 'bottom') +
  theme(legend.position = 'right', legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5))
dev.off()

