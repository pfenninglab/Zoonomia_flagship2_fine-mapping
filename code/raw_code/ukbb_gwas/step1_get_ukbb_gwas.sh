PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping
DATADIR=${PROJDIR}/data/tidy_data
mkdir -p $DATADIR/sumstats $DATADIR/polyfun/munged
cd $DATADIR/sumstats

## download the sumstats
## the tsv table is similar to the xlsx file from https://storage.googleapis.com/broad-alkesgroup-public/UKBB
for ID in {1..47}; do
FILE=$(awk -F '\t' -v VAR=$ID 'FNR==(VAR + 1) {print $2}' ${DATADIR}/tables/readme_ukbb_gwas.tsv )
echo "$ID: $FILE"
if [[ ! -f ${FILE}.sumstats.gz ]]; then wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/${FILE}.sumstats.gz; fi
done

## munge for polyfun format
conda activate polyfun
for ID in {1..34}; do
FILE=$(awk -F '\t' -v VAR=$ID 'FNR==(VAR + 1) {print $2}' ${DATADIR}/tables/readme_ukbb_gwas.tsv )
if [[ $ID -lt 27 ]]; then PREFIX=${FILE}"-Loh_2018"; else PREFIX=${FILE}"-Gazal_2022"; fi
N=$(awk -F '\t' -v VAR=$ID 'FNR==(VAR + 1) {print $4}' ${DATADIR}/tables/readme_ukbb_gwas.tsv )
if [[ ! -f $DATADIR/polyfun/munged/${PREFIX}.parquet ]]; then
# zcat $DATADIR/sumstats/${FILE}.sumstats.gz | \
# awk -v OFS='\t' -v N=$N '{if (NR == 1) {print} else {$11 = N; print}}'| \
# gzip > $DATADIR/sumstats/${FILE}_Neff.sumstats.gz
echo "munging ${FILE}." 
python ~/src/polyfun/munge_polyfun_sumstats.py \
--sumstats $DATADIR/sumstats/${FILE}.sumstats.gz \
--min-info 0.6 --min-maf 0 --n $N \
--out $DATADIR/polyfun/munged/${PREFIX}.parquet
# rm $DATADIR/sumstats/${FILE}_Neff.sumstats.gz
fi
done

