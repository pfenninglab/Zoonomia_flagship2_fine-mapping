SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments

cd $SETDIR

rsync -Pav $GWASDIR/polyfun/munged $SETDIR/data/tidy_data/polyfun
mkdir -p $SETDIR/data/tidy_data/tables
cp $GWASDIR/gwas/unmunged_sumstats/UKBB_BOLTLMM/zreadme_ukbb_gwas.tsv $SETDIR/data/tidy_data/tables/
