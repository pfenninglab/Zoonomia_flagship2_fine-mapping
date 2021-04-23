SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
PHYLOPDIR=${GWASDIR}/phyloP
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
CHAIN=/home/bnphan/resources/liftOver_chainz/hg38ToHg19.over.chain.gz
mkdir -p $DATADIR $DATADIR/bed_hg19 ${SETWD}/data/tidy_data/bed_hg38/HAR $CODEDIR/logs; cd $CODEDIR

######################
## map peaks to hg19
HAR_FILES=$(ls $PHYLOPDIR/phyloP_cutoff_regions/*HARs_202104*.bed.gz )
for FILE in $HAR_FILES; do
NAME=$(basename $FILE | sed 's/.narrowPeak.gz//g;s/.bed.gz//g;s/.hg38//g')
rsync $FILE ${SETWD}/data/tidy_data/bed_hg38/HAR
FILE1500BP=${SETWD}/data/tidy_data/bed_hg38/HAR/$NAME.1500bp.bed.gz
zcat $FILE | awk -F"\t" -v OFS="\t" -v PAD=750 '{MID = ($2 + $3)/2; $2 = MID - PAD; $3 = MID + PAD; print}' | gzip > $FILE1500BP
done

###################################################################
## get the phyloP regions, phastCons regions, C/HARs, ENCODE3 cCREs
ALL_FILES=$(ls ${SETWD}/data/tidy_data/bed_hg38/HAR/*HARs_202104*.bed.gz \
$PHYLOPDIR/phyloP_cutoff_regions/*241mam*.bed.gz \
$PHYLOPDIR/phastCons_cutoff_regions/*43prim*.bed.gz \
${SETWD}/data/tidy_data/bed_hg38/ENCODE3/*.bed.gz | sed '/viterbi/d' )

## map peaks to hg19
for FILE in $ALL_FILES; do
NAME=$(basename $FILE | sed 's/.narrowPeak.gz//g;s/.bed.gz//g;s/.hg38//g')
echo "lifting ${NAME} to hg19"
if [[ ! -f $DATADIR/bed_hg19/${NAME}.hg19.bed.gz ]]; then 
zcat $FILE > $DATADIR/bed_hg19/${NAME}.hg38.bed
liftOver $DATADIR/bed_hg19/${NAME}.hg38.bed $CHAIN $DATADIR/bed_hg19/${NAME}.hg19.bed /dev/null
gzip $DATADIR/bed_hg19/${NAME}.hg19.bed; rm $DATADIR/bed_hg19/${NAME}.hg38.bed
fi
done
