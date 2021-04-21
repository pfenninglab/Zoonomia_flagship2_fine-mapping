SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${GWASDIR}/polyfun/baselineLF2.2.UKB
POLYFUNDIR='/home/bnphan/src/polyfun'
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
ANNOTDIR=${DATADIR}/annotation
cd $CODEDIR; mkdir -p $DATADIR $CODEDIR/logs $DATADIR/annotation

conda activate polyfun

#################################
## 4) annotate the caudate peaks 
LABEL=Stauffer_striatum
for CHR in {1..22}; do

## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
for FILE in $DATADIR/bed/Macaque_striatum*.hg19.bed.gz; do
echo "Extracting chr${CHR} for $(basename $FILE)."
NAME=$(basename $FILE .hg19.bed.gz)
ANNOTFILE2=${ANNOTDIR}/striatum.${NAME}.hg19.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
python $POLYFUNDIR/make_annot_polyfun.py --annot-file ${ANNOTFILE2} --name ${NAME} \
--bed-file ${FILE} --parquet-file ${DATADIR}/UKB_merged_snps.${CHR}.parquet
fi
done

## merge together the annotations together into 1 parquet file per chr
PARQUET=$(ls ${ANNOTDIR}/striatum.*.${CHR}.annot.parquet | sed '/MSN/!d' | tr '\n' ','| sed 's/,$//g')
python $POLYFUNDIR/merge_annot_polyfun.py --annot-file ${ANNOTFILE} \
--parquet-in ${PARQUET},${ANNOTDIR}/caudate_zoonomia_baseline.${CHR}.annot.parquet
fi
done

#####################################################################
## 6) add the base column calculate the M file from the annot parquet
LABEL=Stauffer_striatum
for CHR in {1..22}; do 
python ${POLYFUNDIR}/make_M_polyfun.py --out ${ANNOTDIR}/${LABEL}.${CHR} \
--parquet-file ${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
done

###############################################
## 5) compute LD-scores for phyloP annotations
LABEL=Stauffer_striatum
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=short1 --time 2:00:00 --job-name=${LABEL} --mem=62G --array=14-22 \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"

# rm ${ANNOTDIR}/Stauffer.*
