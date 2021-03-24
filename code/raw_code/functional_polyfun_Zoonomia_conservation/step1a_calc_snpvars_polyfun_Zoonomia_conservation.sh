#!/bin/bash
#SBATCH --partition=pfen_bigmem,gpu,pool3-bigmem,pfen3
#SBATCH --time 12:00:00
#SBATCH --job-name=Asnpvars
#SBATCH --mem=120G
#SBATCH --error=logs/calc_snpvars_%A_%a.txt
#SBATCH --output=logs/calc_snpvars_%A_%a.txt
#SBATCH --array=1-24%1

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/baselineLF2.2.UKB
ZOONOMIADIR=${SETWD}/data/raw_data/zoonomia_annotations/annotation
DATADIR=${SETWD}/data/raw_data/functional_polyfun_Zoonomia_conservation
CODEDIR=${SETWD}/code/raw_code/functional_polyfun_Zoonomia_conservation
POLYFUNDIR='/home/bnphan/src/polyfun'

cd $CODEDIR; source activate polyfun

PREFIX=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
N=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)
CUTOFF=5e-8

echo "Working on ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."

##############################################################################################
# https://github.com/omerwe/polyfun/wiki/1.-Computing-prior-causal-probabilities-with-PolyFun
SUMSTATS=${SETWD}/data/tidy_data/polyfun/munged/${PREFIX}.parquet
OUTDIR=${DATADIR}/${PREFIX}/snpvars
mkdir -p $OUTDIR $DATADIR $CACHEDIR

# 2. Run PolyFun with L2-regularized S-LDSC, using LF2.2.UKB
if [ ! -f ${OUTDIR}/${PREFIX}.22.bins.parquet ]; then
python ${POLYFUNDIR}/polyfun.py --compute-h2-L2 --allow-missing \
--ref-ld-chr ${ZOONOMIADIR}/Zoonomia_conservation. \
--w-ld-chr ${ANNOTDIR}/weights.UKB. \
--output-prefix ${OUTDIR}/$PREFIX --sumstats $SUMSTATS
fi

