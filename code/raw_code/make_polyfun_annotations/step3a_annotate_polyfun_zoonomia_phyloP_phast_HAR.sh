#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
#SBATCH --time=0-8
#SBATCH --mem=20G
#SBATCH --array=1-22
#SBATCH --job-name=zoo_annot
#SBATCH --error=logs/zoo_annot_%A_%a.txt
#SBATCH --output=logs/zoo_annot_%A_%a.txt

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${GWASDIR}/polyfun/baselineLF2.2.UKB
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
POLYFUNDIR='/home/bnphan/src/polyfun'
ANNOTDIR=${DATADIR}/annotation
cd $CODEDIR; mkdir -p $DATADIR $CODEDIR/logs $DATADIR/annotation

source ~/.bashrc; conda activate polyfun

#################################################
## 1) annotate the phyloP scores, phastcons, HARs
LABEL=Zoonomia_annot
FILES=$(ls $DATADIR/bed_hg19/*.hg19.bed.gz | sed '/ENCODE3/d;/HAR/d')
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
for FILE in $FILES ; do
echo "Extracting chr${SLURM_ARRAY_TASK_ID} for $(basename $FILE)."
NAME=$(basename $FILE .hg19.bed.gz | sed 's/^p/Zoonomia_p/g' )
ANNOTFILE2=${ANNOTDIR}/zoo.${NAME}.hg19.${SLURM_ARRAY_TASK_ID}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
TMP=${DATADIR}/bed_hg19/${NAME}.${SLURM_ARRAY_TASK_ID}.tmp.bed.gz
if [[ ! -f $TMP ]]; then zcat $FILE | awk -v SLURM_ARRAY_TASK_ID="chr${SLURM_ARRAY_TASK_ID}" '$1 == SLURM_ARRAY_TASK_ID {print}' | gzip > $TMP; fi
python $POLYFUNDIR/make_annot_polyfun.py --annot-file ${ANNOTFILE2} --name ${NAME} \
--bed-file ${TMP} --parquet-file ${DATADIR}/UKB_merged_snps.${SLURM_ARRAY_TASK_ID}.parquet
rm $TMP
fi; done
## merge together the annotations together into 1 parquet file per chr
PARQUET=$(ls ${ANNOTDIR}/zoo.*.${SLURM_ARRAY_TASK_ID}.annot.parquet | sed '/ENCODE3/d;/HAR/d' | tr '\n' ','| sed 's/,$//g')
python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET} --annot-file ${ANNOTFILE}
fi


#################################################
## 1b) annotate the ENCODE3 annotations separately
FILES=$(ls $DATADIR/bed_hg19/ENCODE3*.hg19.bed.gz )
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet
for FILE in $FILES ; do
echo "Extracting chr${SLURM_ARRAY_TASK_ID} for $(basename $FILE)."
NAME=$(basename $FILE .hg19.bed.gz | sed 's/^p/Zoonomia_p/g' )
ANNOTFILE2=${ANNOTDIR}/encode3.${NAME}.hg19.${SLURM_ARRAY_TASK_ID}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
TMP=${DATADIR}/bed_hg19/${NAME}.${SLURM_ARRAY_TASK_ID}.tmp.bed.gz
if [[ ! -f $TMP ]]; then zcat $FILE | awk -v SLURM_ARRAY_TASK_ID="chr${SLURM_ARRAY_TASK_ID}" '$1 == SLURM_ARRAY_TASK_ID {print}' | gzip > $TMP; fi
python $POLYFUNDIR/make_annot_polyfun.py --annot-file ${ANNOTFILE2} --name ${NAME} \
--bed-file ${TMP} --parquet-file ${DATADIR}/UKB_merged_snps.${SLURM_ARRAY_TASK_ID}.parquet
rm $TMP
fi; done


#####################################################################
## 2) add the base column calculate the M file from the annot parquet
if [[ ! -f ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.l2.M ]]; then python ${POLYFUNDIR}/make_M_polyfun.py --parquet-file ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}.annot.parquet --out ${ANNOTDIR}/${LABEL}.${SLURM_ARRAY_TASK_ID}; fi

###############################################
## 3) compute LD-scores for phyloP annotations
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pfen1,pfen3 --job-name=${LABEL} --mem=45G --array=${SLURM_ARRAY_TASK_ID} --time 1-0:00:00 \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"

# rm ${ANNOTDIR}/zoo.* ${ANNOTDIR}/encode3.*



