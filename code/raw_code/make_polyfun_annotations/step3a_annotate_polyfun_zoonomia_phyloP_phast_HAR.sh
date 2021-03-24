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
## 1) annotate the phyloP scores, phastcons, HARs
LABEL=Zoonomia_phyloP_phast_HAR
FILES=$(ls $DATADIR/bed/phyloPam*.hg19.bed.gz $DATADIR/bed/PhastCons*.hg19.bed.gz $DATADIR/bed/zoonomia_HARs_20210304.hg19.bed.gz | sed '/extend/d')
for CHR in {1..22}; do
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
for FILE in $FILES ; do
echo "Extracting chr${CHR} for $(basename $FILE)."
NAME=Zoonomia_$(basename $FILE .hg19.bed.gz | sed 's/zoonomia_HARs_20210304/HAR/g;s/_viterbi_conservedRegions//g')
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


#####################################################################
## 2) annotate the 500bp extensions phyloP scores, phastcons, HARs
LABEL=Zoonomia_phyloP_phast_HAR_extend500
FILES=$(ls $DATADIR/bed/phyloPam*.extend500.hg19.bed.gz $DATADIR/bed/PhastCons*.extend500.hg19.bed.gz $DATADIR/bed/zoonomia_HARs_20210304.extend500.hg19.bed.gz)
for CHR in {1..22}; do
## the final annotation file
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then
for FILE in $FILES ; do
echo "Extracting chr${CHR} for $(basename $FILE)."
NAME=Zoonomia_$(basename $FILE .hg19.bed.gz | sed 's/zoonomia_HARs_20210304/HAR/g;s/_viterbi_conservedRegions//g')
ANNOTFILE2=${ANNOTDIR}/ext500.${NAME}.hg19.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE2 ]]; then
TMP=${DATADIR}/bed/${NAME}.${CHR}.tmp.bed.gz
if [[ ! -f $TMP ]]; then zcat $FILE | awk -v CHR="chr${CHR}" '$1 == CHR {print}' | gzip > $TMP; fi
python $POLYFUNDIR/make_annot_polyfun.py --annot-file ${ANNOTFILE2} --name ${NAME} \
--bed-file ${TMP} --parquet-file ${DATADIR}/UKB_merged_snps.${CHR}.parquet
rm $TMP
fi; done

## merge together the annotations together into 1 parquet file per chr
ANNOTFILE1=${ANNOTDIR}/Zoonomia_phyloP_extend500.${CHR}.annot.parquet
PARQUET1=$(ls ${ANNOTDIR}/ext500.*phyloPam*.${CHR}.annot.parquet | tr '\n' ','| sed 's/,$//g')
if [[ ! -f $ANNOTFILE1 ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET1} --annot-file ${ANNOTFILE1}; fi

ANNOTFILE2=${ANNOTDIR}/Zoonomia_phast_HAR_extend500.${CHR}.annot.parquet
PARQUET2=$(ls ${ANNOTDIR}/ext500.*PhastCons*.${CHR}.annot.parquet ${ANNOTDIR}/ext500.*HAR*.${CHR}.annot.parquet | tr '\n' ','| sed 's/,$//g')
if [[ ! -f $ANNOTFILE2 ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET2} --annot-file ${ANNOTFILE2}; fi

## merge together the annotations together into 1 parquet file per chr
if [[ -f $ANNOTFILE1 &&  -f $ANNOTFILE2 ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${ANNOTFILE1},${ANNOTFILE2} --annot-file ${ANNOTFILE}; fi
fi; done


############################################################
# 3) merge the orig annotations together w/ the extended ###
LABEL=Zoonomia_conservation
for CHR in {1..10}; do
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
## merge together the annotations together into 1 parquet file per chr
if [[ ! -f $ANNOTFILE ]]; then 
PARQUET=${ANNOTDIR}/Zoonomia_phyloP_phast_HAR.${CHR}.annot.parquet,${ANNOTDIR}/Zoonomia_phyloP_phast_HAR_extend500.${CHR}.annot.parquet
python $POLYFUNDIR/merge_annot_polyfun.py \
--parquet-in ${PARQUET} --annot-file ${ANNOTFILE}; 
fi; done


#####################################################################
## 3) add the base column calculate the M file from the annot parquet
LABEL=Zoonomia_conservation
for CHR in {1..10}; do 
if [[ ! -f ${ANNOTDIR}/${LABEL}.${CHR}.l2.M ]]; then python ${POLYFUNDIR}/make_M_polyfun.py --parquet-file ${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet --out ${ANNOTDIR}/${LABEL}.${CHR}; fi
done

###############################################
## 2) compute LD-scores for phyloP annotations
LABEL=Zoonomia_conservation
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pfen3 --time 23:00:00 --job-name=${LABEL} --mem=60G --array=1-10 --time 1-0:00:00 \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"



