SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CACHEDIR=${GWASDIR}/polyfun/LD_cache
ANNOTDIR=${GWASDIR}/polyfun/baselineLF2.2.UKB
PHYLOPDIR=${GWASDIR}/phyloP/phyloP_cutoff_regions
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
CODEDIR=${SETWD}/code/raw_code/make_polyfun_annotations
CHAIN=/home/bnphan/resources/liftOver_chainz/hg38ToHg19.over.chain.gz
CAUDATE_PEAKS=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/hg38/Corces_2020/peak
MONKEY_PEAKS=/projects/pfenninggroup/singleCell/Macaque_snATAC-seq/macaque_snATAC-seq/m015_Peanut_snATAC/data/raw_data/halper
mkdir -p $DATADIR $DATADIR/bed $CODEDIR/logs; cd $CODEDIR

#############################################################
## get the phyloP scores, Corces2020 caudate scATAC peakset
ALL_FILES=$(ls $CAUDATE_PEAKS/*narrowPeak.gz $PHYLOPDIR/*.bed.gz \
	$MONKEY_PEAKS/*Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz | \
	sed '/Consensus/d;/mappedTo/d;/UNK/d' )

######################
## map peaks to hg19
for FILE in $ALL_FILES; do
NAME=$(basename $FILE | sed 's/.narrowPeak.gz//g;s/.bed.gz//g;s/.GenBankRheMac8.Macaca_mulattaToHomo_sapiens.HALPER//g;s/Stauffer_striatum/Macaque_striatum/g;s/.hg38//g')
if [[ ! -f $DATADIR/bed/${NAME}.hg19.bed.gz ]]; then 
zcat $FILE > $DATADIR/bed/${NAME}.hg38.bed
liftOver $DATADIR/bed/${NAME}.hg38.bed $CHAIN $DATADIR/bed/${NAME}.hg19.bed /dev/null
gzip $DATADIR/bed/${NAME}.hg19.bed; rm $DATADIR/bed/${NAME}.hg38.bed
fi
done

##############################################################
## extend Zoonomia annotation peaks on either sides by 500bp
HG19_CHR=/home/bnphan/resources/genomes/hg19/hg19.chrom.sizes
FILES=$(ls $DATADIR/bed/*.hg19.bed.gz | sed '/extend/d;/Corces2020/d;/Macaque_striatum/d' )
for FILE in ; do
EXTFILE=$DATADIR/bed/$(basename $FILE .hg19.bed.gz).extend500.hg19.bed
if [[ ! -f ${EXTFILE}.gz ]]; then 
TMP=$DATADIR/bed/$(basename $FILE .gz).tmp
TMP2=$DATADIR/bed/$(basename $FILE .gz).tmp2
TMP3=$DATADIR/bed/$(basename $FILE .gz).tmp3
zcat $FILE > $TMP
bedtools slop  -i $TMP -g $HG19_CHR -b 500 > $TMP2
sort -k1,1 -k2,2n $TMP2 > $TMP3
bedtools merge -i $TMP3 > $EXTFILE
gzip $EXTFILE; rm $TMP $TMP2 $TMP3
fi; done


######################################################
### work on the zoonomia human accelerated regions ###
## zoomHARs_named.bed received from Kathleen Keough March 3rd, 2021
## liftover hg38 coordinates to hg19
# cat $DATADIR/bed/zoomHARs_named.bed | sort -k 1,1 -k 2,2n | gzip > $DATADIR/bed/zoonomia_HARs_20210304.bed.gz
# liftOver $DATADIR/bed/zoomHARs_named.bed $CHAIN $DATADIR/bed/zoonomia_HARs_20210304.hg19.bed /dev/null
# gzip $DATADIR/bed/zoonomia_HARs_20210304.hg19.bed
