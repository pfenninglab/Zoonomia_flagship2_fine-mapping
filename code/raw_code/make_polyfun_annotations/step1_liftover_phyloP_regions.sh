SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CACHEDIR=${GWASDIR}/polyfun/LD_cache
ANNOTDIR=${GWASDIR}/polyfun/baselineLF2.2.UKB
PHYLOPDIR=${GWASDIR}/phyloP/phyloP_cutoff_regions
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
CHAIN=/home/bnphan/resources/liftOver_chainz/hg38ToHg19.over.chain.gz
CAUDATE_PEAKS=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/hg38/Corces_2020/peak
mkdir -p $DATADIR $DATADIR/bed $CODEDIR/logs; cd $CODEDIR

## get the phyloP scores, Corces2020 caudate scATAC peakset
ALL_FILES=$(ls $CAUDATE_PEAKS/*narrowPeak.gz \
	$PHYLOPDIR/200m_scoresPhyloP_20210214.*.05.bed.gz \
	$PHYLOPDIR/200m_scoresPhyloP_20210214.*.05.bed.gz | \
	sed '/Consensus/d;/mappedTo/d')

for FILE in $ALL_FILES; do
NAME=$(basename $FILE | sed 's/.narrowPeak.gz//g;s/.bed.gz//g')
if [[ ! -f $DATADIR/bed/${NAME}.hg19.bed.gz ]]; then 
zcat $FILE > $DATADIR/bed/${NAME}.hg38.bed
liftOver $DATADIR/bed/${NAME}.hg38.bed $CHAIN $DATADIR/bed/${NAME}.hg19.bed /dev/null
gzip $DATADIR/bed/${NAME}.hg19.bed; rm $DATADIR/bed/${NAME}.hg38.bed
fi
done

######################################################
### work on the zoonomia human accelerated regions ###
## zoomHARs_named.bed received from Kathleen Keough March 3rd, 2021
## liftover hg38 coordinates to hg19
cat $DATADIR/bed/zoomHARs_named.bed | sort -k 1,1 -k 2,2n | gzip > $DATADIR/bed/zoonomia_HARs_20210304.bed.gz
liftOver $DATADIR/bed/zoomHARs_named.bed $CHAIN $DATADIR/bed/zoonomia_HARs_20210304.hg19.bed /dev/null
gzip $DATADIR/bed/zoonomia_HARs_20210304.hg19.bed





