#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1,pfen3,pfen_bigmem
#SBATCH --time=0-8
#SBATCH --mem=20G
#SBATCH --array=1-22
#SBATCH --job-name=merge_annot
#SBATCH --error=logs/merge_annot_%A_%a.txt
#SBATCH --output=logs/merge_annot_%A_%a.txt

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
POLYFUNDIR='/home/bnphan/src/polyfun'
ANNOTDIR=${DATADIR}/annotation
cd $CODEDIR; mkdir -p $DATADIR $CODEDIR/logs $DATADIR/annotation

source ~/.bashrc; conda activate polyfun

#####################################
## 1) merge the annotations from the 
LABEL=Zoonomia_annot_baselineLF
ANNOTFILE=${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet
PARQUET=${ANNOTDIR}/Zoonomia_annot.${SLURM_ARRAY_TASK_ID}.annot.parquet,${GWASDIR}/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.${SLURM_ARRAY_TASK_ID}.annot.parquet
PARQUET2=$(ls ${ANNOTDIR}/encode3*.${SLURM_ARRAY_TASK_ID}.annot.parquet| sed '/HAR/d' | tr '\n' ','| sed 's/,$//g')
if [[ ! -f $ANNOTFILE ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET},${PARQUET2} --annot-file ${ANNOTFILE} --drop-regex 'Conserved_LindbladToh|Conserved_Mammal|Conserved_Primate'
fi

######################################################################
## 2) add the base column calculate the M file from the annot parquet
if [[ ! -f ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.l2.M ]]; then python ${POLYFUNDIR}/make_M_polyfun.py --parquet-file ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet --out ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}; fi

###############################################
## 3) compute LD-scores for phyloP annotations
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pfen1 --time 1-00:00:00 --job-name=merge --mem=45G --array=$SLURM_ARRAY_TASK_ID \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"

# rm ${ANNOTDIR}/zoo.*.${SLURM_ARRAY_TASK_ID}.annot.parquet
