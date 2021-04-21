SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
POLYFUNDIR='/home/bnphan/src/polyfun'
ANNOTDIR=${DATADIR}/annotation
cd $CODEDIR; mkdir -p $DATADIR $CODEDIR/logs $DATADIR/annotation

conda activate polyfun


#####################################
## 1) merge the annotations from the 
LABEL=Zoonomia_conservation_baselineLF
for CHR in {1..22}; do 
ANNOTFILE=${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet
PARQUET=${ANNOTDIR}/Zoonomia_conservation.${CHR}.annot.parquet,${GWASDIR}/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.${CHR}.annot.parquet
if [[ ! -f $ANNOTFILE ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --parquet-in ${PARQUET} --annot-file ${ANNOTFILE} --drop-regex 'Conserved_LindbladToh|Conserved_Mammal|Conserved_Primate'
fi; done


######################################################################
## 2) add the base column calculate the M file from the annot parquet
LABEL=Zoonomia_conservation_baselineLF
for CHR in {1..22}; do 
if [[ ! -f ${ANNOTDIR}/${LABEL}.${CHR}.l2.M ]]; then python ${POLYFUNDIR}/make_M_polyfun.py --parquet-file ${ANNOTDIR}/${LABEL}.${CHR}.annot.parquet --out ${ANNOTDIR}/${LABEL}.${CHR}; fi
done


###############################################
## 3) compute LD-scores for phyloP annotations
LABEL=Zoonomia_conservation_baselineLF
THECALL="source ~/.bashrc; conda activate polyfun; if [[ ! -f ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ]]; then python ${POLYFUNDIR}/compute_ldscores_ukb.py --annot ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.annot.parquet --ld-dir ${GWASDIR}/polyfun/LD_cache --out ${ANNOTDIR}/${LABEL}.\${SLURM_ARRAY_TASK_ID}.l2.ldscore.parquet ; fi"
sbatch --partition=pfen1 --time 2:00:00 --job-name=merge --mem=45G --array=1-22 \
--error=logs/run_annot_%A_%a.out.txt --output=logs/run_annot_%A_%a.out.txt --wrap="${THECALL}"

