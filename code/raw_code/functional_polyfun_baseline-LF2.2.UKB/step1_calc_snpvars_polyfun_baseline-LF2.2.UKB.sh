#!/bin/bash
#SBATCH --partition=pfen_bigmem
#SBATCH --time 4-4
#SBATCH --job-name=snpvars
#SBATCH --mem=60G
#SBATCH --error=logs/calc_snpvars_%A_%a.txt
#SBATCH --output=logs/calc_snpvars_%A_%a.txt
#SBATCH --array=24

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
ANNOTDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/baselineLF2.2.UKB
DATADIR=${SETWD}/data/raw_data/functional_polyfun_baseline-LF2.2.UKB
CODEDIR=${SETWD}/code/raw_code/functional_polyfun_baseline-LF2.2.UKB
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
--ref-ld-chr ${ANNOTDIR}/baselineLF2.2.UKB. \
--w-ld-chr ${ANNOTDIR}/weights.UKB. \
--output-prefix ${OUTDIR}/$PREFIX --sumstats $SUMSTATS
fi

# 3. Compute LD-scores for each SNP bin, do per chromosome
for CHR in {1..22}; do
if [ ! -f ${OUTDIR}/${PREFIX}.${CHR}.l2.ldscore.parquet ]; then
python ${POLYFUNDIR}/polyfun.py --compute-ldscores --chr ${CHR} \
--output-prefix ${OUTDIR}/$PREFIX --ld-ukb --ld-dir $CACHEDIR
fi
done

# 4. Re-estimate per-SNP heritabilities via S-LDSC
if [ ! -f ${OUTDIR}/${PREFIX}.22.snpvar_constrained.gz ]; then
python ${POLYFUNDIR}/polyfun.py --compute-h2-bins \
--allow-missing --output-prefix ${OUTDIR}/$PREFIX \
--sumstats $SUMSTATS --w-ld-chr ${ANNOTDIR}/weights.UKB.
fi



