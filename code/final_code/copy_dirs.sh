SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments

cd $SETDIR

# rsync -Pav $GWASDIR/polyfun/munged $SETDIR/data/tidy_data/polyfun
# mkdir -p $SETDIR/data/tidy_data/tables
# cp $GWASDIR/gwas/unmunged_sumstats/UKBB_BOLTLMM/zreadme_ukbb_gwas.tsv $SETDIR/data/tidy_data/tables/


## upload to BU_Addiction_snRNA-seq/Striatum
https://drive.google.com/drive/folders/1aep9c9bIPAf3itNsz2apOnA6T0QbPFv6?usp=sharing
DATADIR=$SETDIR/data/raw_data/zoonomia_annotations/annotation
gdrive upload --recursive --parent 1aep9c9bIPAf3itNsz2apOnA6T0QbPFv6 $DATADIR/to_upload
