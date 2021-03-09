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
## 1) annotate the phyloP scores
LABEL=200m_scoresPhyloP_20210214_HAR_20210304
for CHR in {1..8}; do
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
for FILE in $DATADIR/bed/200m_scoresPhyloP_20210214.*.hg19.bed.gz $DATADIR/bed/zoonomia_HARs_20210304.hg19.bed.gz; do
echo "Extracting chr${CHR} for $(basename $FILE)."
NAME=$(basename $FILE .hg19.bed.gz | sed 's/200m_scoresPhyloP_20210214./phyloP./g;s/zoonomia_HARs_20210304/zoomHAR/g')
ANNOTFILE2=${ANNOTDIR}/241m.${NAME}.hg19.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
TMP=${DATADIR}/bed/${NAME}.${CHR}.tmp.bed.gz
if [[ ! -f $TMP ]]; then zcat $FILE | awk -v CHR="chr${CHR}" '$1 == CHR {print}' | gzip > $TMP; fi
python $POLYFUNDIR/make_annot_polyfun.py --annot-file ${ANNOTFILE2} --name ${NAME} \
--bed-file ${TMP} --parquet-file ${DATADIR}/UKB_merged_snps.${CHR}.parquet
rm $TMP
fi; done
## merge together the annotations together into 1 parquet file per chr
PARQUET=$(ls ${ANNOTDIR}/241m.*.${CHR}.annot.parquet | tr '\n' ','| sed 's/,$//g')
python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET} --annot-file ${ANNOTFILE}
fi; done

rm ${ANNOTDIR}/241m.*


###############################################
## 2) compute LD-scores for phyloP annotations
LABEL=200m_scoresPhyloP_20210214_HAR_20210304
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pool1 --time 24:00:00 --job-name=${LABEL} --mem=43G --array=8,16-22 --time 1-0:00:00 \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"


#####################################################################
## 3) add the base column calculate the M file from the annot parquet
LABEL=200m_scoresPhyloP_20210214_HAR_20210304
for CHR in {1..22}; do 
if [[ ! -f ${ANNOTDIR}/${LABEL}.${CHR}.l2.M ]]; then python ${POLYFUNDIR}/make_M_polyfun.py --parquet-file ${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet --out ${ANNOTDIR}/${LABEL}.${CHR}; fi
done

