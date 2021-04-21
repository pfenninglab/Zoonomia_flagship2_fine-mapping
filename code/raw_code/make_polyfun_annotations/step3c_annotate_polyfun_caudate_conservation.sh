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
LABEL=caudate_zoonomia_baseline
for CHR in {1..22}; do
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
for FILE in $DATADIR/bed/Corces2020_caudate*.hg19.bed.gz; do
echo "Extracting chr${CHR} for $(basename $FILE)."
NAME=$(basename $FILE .hg19.bed.gz | sed 's/Corces2020_caudate/Caud/g')
ANNOTFILE2=${ANNOTDIR}/${NAME}.hg19.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
python $POLYFUNDIR/make_annot_polyfun.py --annot-file ${ANNOTFILE2} --name ${NAME} \
--bed-file ${FILE} --parquet-file ${DATADIR}/UKB_merged_snps.${CHR}.parquet
fi
done
## merge together the annotations together into 1 parquet file per chr
PARQUET2=$(ls ${ANNOTDIR}/Caud.*.${CHR}.annot.parquet | tr '\n' ','| sed 's/,$//g')
ANNOTFILE2=${ANNOTDIR}/Caud.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET2} --annot-file ${ANNOTFILE2}; fi
## merge together the the rhesus-human orth annotations together into 1 parquet file per chr
PARQUET3=$(ls ${ANNOTDIR}/Caud_hgRmOrth.*.${CHR}.annot.parquet | tr '\n' ','| sed 's/,$//g')
ANNOTFILE3=${ANNOTDIR}/Caud_hgRmOrth.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE3 ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET3} --annot-file ${ANNOTFILE3}; fi
## merge together the mouse-human orth annotations together into 1 parquet file per chr
PARQUET4=$(ls ${ANNOTDIR}/Caud_hgMmOrth.*.${CHR}.annot.parquet | tr '\n' ','| sed 's/,$//g')
ANNOTFILE4=${ANNOTDIR}/Caud_hgMmOrth.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE4 ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET4} --annot-file ${ANNOTFILE4}; fi
# merge all of the annotations together
python $POLYFUNDIR/merge_annot_polyfun.py --annot-file ${ANNOTFILE} \
--parquet-in ${ANNOTFILE2},${ANNOTFILE3},${ANNOTFILE4},${ANNOTDIR}/Zoonomia_conservation_baselineLF.${CHR}.annot.parquet
fi
done


#####################################################################
## 6) add the base column calculate the M file from the annot parquet
LABEL=caudate_zoonomia_baseline
for CHR in {1..22}; do 
python ${POLYFUNDIR}/make_M_polyfun.py --out ${ANNOTDIR}/${LABEL}.${CHR} \
--parquet-file ${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
done


###############################################
## 5) compute LD-scores for phyloP annotations
LABEL=caudate_zoonomia_baseline
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pool3-bigmem --time 12:00:00 --job-name=${LABEL} --mem=62G --array=1-22 \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"

# rm ${ANNOTDIR}/Caud*
