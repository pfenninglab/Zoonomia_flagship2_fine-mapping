# Fine-mapping human complex traits with mammalian genome constraint
Code to run polyfun functional fine-mapping with Zoonomia annotations.

## Table of contents:
The following code directories are organized by various taskes

code/raw_code:
- ukbb_gwas: download and pre-process the Loh et al. 2018 GWAS
- nonfunct_finemapping: performn non-functional fine-mapping with SuSIE
- make_polyfun_annotations: annotate PolyFun baseline and add Zoonomia annotations
- functional_polyfun_Zoonomia_phyloP_phast_HAR: run PolyFun with only Zoonomia annotations
- functional_polyfun_baseline-LF2.2.UKB: run PolyFun with the BaselineLF2 annotations
- functional_polyfun_Zoonomia_annot_baselineLF2: run PolyFun with the combined annotations

figures/explanatory:
- finemapping_figures: cumulative distribution function plots for Zoonomia Flagship 2 manuscript

figures/exploratory:
-comparing_improved_PIP_scores: compare fine-mapped SNPs w/ UNICORNS
-functional_polyfun_Zoonomia_phyloP: compare fine-mapped SNPs across annotations

## Data:
The following describes the outputs from PolyFun aggregating fine-mapped SNPs across traits
data/raw_data/functional_polyfun_Zoonomia_annot_baselineLF2:
- "\_top_annot.txt.gz": fine-mapped SNPs along with corresponding annotations
- "\_causal_set.txt.gz": fine-mapped SNPs with non-zero PIP
