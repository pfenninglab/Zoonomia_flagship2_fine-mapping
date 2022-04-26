#!/bin/bash
#SBATCH --partition=short1,gpu,pool1
#SBATCH --job-name=ukb_nonfunct
#SBATCH --mem=20G
#SBATCH --error=logs/polyfun_nonfunct_%A_%a.txt
#SBATCH --output=logs/polyfun_nonfunct_%A_%a.txt
#SBATCH --array=1-34

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
DATADIR=${SETWD}/data/raw_data/nonfunct_finemapping
CODEDIR=${SETWD}/code/raw_code/nonfunct_finemapping
POLYFUNDIR='/home/bnphan/src/polyfun'

cd $CODEDIR; 
source activate polyfun

if [[ $SLURM_ARRAY_TASK_ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
N=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)
CUTOFF=5e-8

echo "Working on ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."

################################################################
# https://github.com/omerwe/polyfun/wiki/3.-Functionally-informed-fine-mapping-with-finemapper
# generate scripts to run poly fun mapping at significant loci
SUMSTATS=${SETWD}/data/tidy_data/polyfun/munged/${PREFIX}.parquet
OUTDIR=${DATADIR}/${PREFIX}/susie
mkdir -p $OUTDIR $DATADIR $CACHEDIR

for CHR in {1..22}; do
## if jobs not already created
if [ ! -f ${OUTDIR}/polyfun_all_jobs_${CHR}.txt ]; then
echo "Working on chr${CHR}."
python ${POLYFUNDIR}/create_finemapper_jobs.py \
--non-funct --allow-missing --method susie \
--cache-dir ${CACHEDIR} --sumstats ${SUMSTATS} \
--pvalue-cutoff ${CUTOFF} --max-num-causal 10 \
--out-prefix ${OUTDIR}/polyfun_all --n ${N} \
--jobs-file ${OUTDIR}/polyfun_all_jobs_${CHR}.txt \
--memory 15 --chr ${CHR}
fi
done
