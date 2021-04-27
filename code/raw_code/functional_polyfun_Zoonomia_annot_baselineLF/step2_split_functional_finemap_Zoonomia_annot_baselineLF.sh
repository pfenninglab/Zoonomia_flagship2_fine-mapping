#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --time 2:00:00
#SBATCH --job-name=split_jobs
#SBATCH --mem=10G
#SBATCH --error=logs/polyfun_funct_Zoonomia_annot_baselineLF_%A_%a.txt
#SBATCH --output=logs/polyfun_funct_Zoonomia_annot_baselineLF_%A_%a.txt
#SBATCH --array=1-24

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
DATADIR=${SETWD}/data/raw_data/functional_polyfun_Zoonomia_annot_baselineLF
CODEDIR=${SETWD}/code/raw_code/functional_polyfun_Zoonomia_annot_baselineLF
POLYFUNDIR='/home/bnphan/src/polyfun'

cd $CODEDIR; source ~/.bashrc; conda activate polyfun

# for SLURM_ARRAY_TASK_ID in {12..15}; do
PREFIX=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
N=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)
CUTOFF=5e-8

echo "Working on ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."

################################################################
# https://github.com/omerwe/polyfun/wiki/3.-Functionally-informed-fine-mapping-with-finemapper
# generate scripts to run poly fun mapping at significant loci
OUTDIR=${DATADIR}/${PREFIX}/susie
mkdir -p $OUTDIR $DATADIR $CACHEDIR

for CHR in {1..22}; do
## if jobs not already created
SUMSTATS=${DATADIR}/${PREFIX}/snpvars/${PREFIX}.${CHR}.snpvar_constrained.gz
if [[ ! -f ${OUTDIR}/polyfun_all_jobs_${CHR}.txt && -f ${DATADIR}/${PREFIX}/snpvars/${PREFIX}.${CHR}.snpvar_constrained.gz ]]; then
echo "Working on chr${CHR}."
python ${POLYFUNDIR}/create_finemapper_jobs.py \
--allow-missing --method susie \
--cache-dir ${CACHEDIR} --sumstats ${SUMSTATS} \
--pvalue-cutoff ${CUTOFF} --max-num-causal 10 \
--out-prefix ${OUTDIR}/polyfun_all --n ${N} \
--jobs-file ${OUTDIR}/polyfun_all_jobs_${CHR}.txt \
--memory 5 --chr ${CHR}
fi
done
# done
