#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1,pfen_bigmem
#SBATCH --time=0-8
#SBATCH --mem=20G
#SBATCH --array=1-22
#SBATCH --job-name=caud_merge_annot
#SBATCH --error=logs/caud_merge_annot_%A_%a.txt
#SBATCH --output=logs/caud_merge_annot_%A_%a.txt

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${GWASDIR}/polyfun/baselineLF2.2.UKB
POLYFUNDIR='/home/bnphan/src/polyfun'
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
ANNOTDIR=${DATADIR}/annotation
cd $CODEDIR; mkdir -p $DATADIR $CODEDIR/logs $DATADIR/annotation

source ~/.bashrc; conda activate polyfun

#################################
## 1) annotate the caudate peaks 
LABEL=Caudate_Zoonomia_annot_baselineLF
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
echo "Dropping old Zoonomia annotations chr${SLURM_ARRAY_TASK_ID} with regex: Zoonomia."
ANNOTFILE2=${ANNOTDIR}/tmpCaudate.hg19.${SLURM_ARRAY_TASK_ID}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
python $POLYFUNDIR/merge_annot_polyfun.py --annot-file ${ANNOTFILE2} \
--parquet-in ${ANNOTDIR}/caudate_zoonomia_baseline.${SLURM_ARRAY_TASK_ID}.annot.parquet \
--drop-regex 'Zoonomia|HAR'
fi

# merge together caudate annotations and new Zoonomia annotations, and ENCODE3 annotations
ANNOTFILE3=$(ls ${ANNOTDIR}/encode3*.${SLURM_ARRAY_TASK_ID}.annot.parquet | tr '\n' ','| sed 's/,$//g')
ANNOTFILE4=${ANNOTDIR}/Zoonomia_annot.${SLURM_ARRAY_TASK_ID}.annot.parquet
python $POLYFUNDIR/merge_annot_polyfun.py --annot-file ${ANNOTFILE} \
--parquet-in ${ANNOTFILE2},${ANNOTFILE3},${ANNOTFILE4}
fi

#####################################################################
## 6) add the base column calculate the M file from the annot parquet
python ${POLYFUNDIR}/make_M_polyfun.py --out ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID} \
--parquet-file ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet

###############################################
## 5) compute LD-scores for phyloP annotations
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pfen1,pfen_bigmem --time 12:00:00 --job-name=${LABEL} --mem=62G --array=${SLURM_ARRAY_TASK_ID} \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"