#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --job-name=aggregate
#SBATCH --mem=30G
#SBATCH --error=logs/runAggregate_%A.txt
#SBATCH --output=logs/runAggregate_%A.txt

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
ANNOTDIR=${DATADIR}/annotation
POLYFUNDIR='/home/bnphan/src/polyfun'
PIP_CUTOFF=0.95
source activate polyfun

###########################################################
# aggregate finemapping results over different SNP blocks
FINESNP=${OUTDIR}/../${PREFIX}_aggregate.txt.gz
if [[ ! -f $FINESNP ]]; then #
python ${POLYFUNDIR}/aggregate_finemapper_results.py --allow-missing-jobs \
--out-prefix ${OUTDIR}/polyfun_all --sumstats ${SUMSTATS} \
--pvalue-cutoff ${CUTOFF} --out ${FINESNP}
fi

################################################
# filter out SNPs that are in the credible set
PIPSNPS=${OUTDIR}/../${PREFIX}_causal_set.txt.gz
echo "Extracting top PIP SNPs on ${CHR} with PIP greater than ${PIP_CUTOFF}."
if [[ ! -f $PIPSNPS && -f $FINESNP ]]; then
# credible set is in last column
zcat $FINESNP | awk -v VAR=${PIP_CUTOFF} '{if( $(NF - 3) >= VAR ) print }' | gzip > $PIPSNPS
fi

##############################################
# filter and annotate top PIP snps per chunks
ANONTSNP=${OUTDIR}/../${PREFIX}_top_annot.txt.gz
if [[ ! -f $ANONTSNP && -f $FINESNP ]]; then
for CHR in {1..22}; do
echo "Annotating SNPs on ${CHR} with PIP greater than ${PIP_CUTOFF}."
FILE=${OUTDIR}/polyfun_all.chr${CHR}.gz
OUT=$(echo $FILE | sed 's/polyfun_all/top_annot/g')
if [[ ! -f $OUT ]]; then 
# find which column is CHR
COL=$(zcat $FINESNP | awk -v col=CHR 'NR==1{for(i=1;i<=NF;i++) {if($i==col) {print i ;exit}} }')
zcat $FINESNP | awk  -v COL=$COL -v VAR=${CHR} '{if( NR == 1 || $COL == VAR ) print }' | gzip > $FILE
python ${POLYFUNDIR}/extract_annotations.py \
--annot ${ANNOTDIR}/merged_baselineLF_200m_conservation.${CHR}.annot.parquet \
--pips $FILE --pip-cutoff $PIP_CUTOFF --out $OUT
rm $FILE
fi
done
cat ${OUTDIR}/top_annot.*.gz | zcat | awk '{if( NR == 1 || $2 ~ /^[0-9]+$/ ) print }' | gzip > $ANONTSNP
rm ${OUTDIR}/top_annot.*.gz
fi
