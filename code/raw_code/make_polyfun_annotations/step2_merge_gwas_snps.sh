SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${GWASDIR}/polyfun/baselineLF2.2.UKB
POLYFUNDIR='/home/bnphan/src/polyfun'
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
ANNOTDIR=${DATADIR}/annotation
cd $CODEDIR; mkdir -p $DATADIR $CODEDIR/logs $DATADIR/annotation

conda activate polyfun

#########################################
## 1) merge the SNPs across parquet files
PARQUET=$(ls ${GWASDIR}/polyfun/baselineLF2.2.UKB/weights.UKB.*.parquet $SETWD/data/tidy_data/polyfun/munged/*.parquet /home/bnphan/projects/Addiction_MPRA_2021/data/raw_data/polyfun/munged/*.parquet | tr '\n' ',' | sed 's/,$//g')
ANNOTFILE=${DATADIR}/UKB_merged_snps.parquet
if [[ ! -f ${ANNOTFILE} ]]; then python $POLYFUNDIR/merge_annot_polyfun.py --snp-only --parquet-in ${PARQUET} --annot-file ${ANNOTFILE}; fi

## 2) split the SNPs by chromosome
if [[ ! -f ${DATADIR}/UKB_merged_snps.22.parquet ]]; then python $POLYFUNDIR/make_M_polyfun.py --split-chr --parquet-file ${DATADIR}/UKB_merged_snps.parquet --out ${DATADIR}/UKB_merged_snps ; fi
