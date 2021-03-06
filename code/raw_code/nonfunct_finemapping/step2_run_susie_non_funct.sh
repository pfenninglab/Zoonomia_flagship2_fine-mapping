SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
DATADIR=${SETWD}/data/raw_data/nonfunct_finemapping
CODEDIR=${SETWD}/code/raw_code/nonfunct_finemapping
POLYFUNDIR='/home/bnphan/src/polyfun'
CUTOFF=5e-8; cd $CODEDIR;


######################################
## split finemapping job by chromosome
for ID in {1..34}; do
if [[ $ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
N=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)
echo "Working on ${PREFIX} GWAS with P -lt ${CUTOFF} cutoff for loci."
SUMSTATS=${SETWD}/data/tidy_data/polyfun/munged/${PREFIX}.parquet
OUTDIR=${DATADIR}/${PREFIX}/susie
mkdir -p $OUTDIR $DATADIR $CACHEDIR
## use wildcard character @ character
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_nonfunct_finemapping_aggregate.txt.gz ]]; then
sbatch -p pfen1 --mem 90G --time 3-00:00:00 --array 1-22 \
--export=JobsFileName="${OUTDIR}/polyfun_all_jobs_@.txt" \
${CODEDIR}/slurm_finemap_byLine.sh
fi
done


###########################################
## combine all the jobs together per trait
## catch the jobs that need a lot of RAM
for ID in {1..34}; do
if [[ $ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
N=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)
OUTDIR=${DATADIR}/${PREFIX}/susie
JobsFileName=${OUTDIR}/polyfun_all_jobs.txt
if [[ ! -f $JobsFileName || $(wc -l $JobsFileName | cut -d ' ' -f1 ) == 0  ]]; then cat ${OUTDIR}/polyfun_all_jobs_*.txt > $JobsFileName; fi
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_nonfunct_finemapping_aggregate.txt.gz ]]; then
echo "Working on ${PREFIX} GWAS with P -lt ${CUTOFF} cutoff for loci."
sbatch -p pool3-bigmem,pfen_bigmem,pfen1 --mem 220G --time 4-04:00:00 --export=JobsFileName=${JobsFileName} ${CODEDIR}/slurm_finemap_byLine.sh
fi
done


##################################################
## submit job to aggregate fine-mapping when done
for ID in {1..34}; do
if [[ $ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
N=$( awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv )
CUTOFF=5e-8
SUMSTATS=${SETWD}/data/tidy_data/polyfun/munged/${PREFIX}.parquet
OUTDIR=${DATADIR}/${PREFIX}/susie
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_nonfunct_finemapping_top_annot.txt.gz ]]; then
echo "Aggregating results from ${PREFIX} GWAS with P -lt ${CUTOFF} cutoff for loci."
sbatch --export=OUTDIR=${OUTDIR},PREFIX=${PREFIX}_nonfunct_finemapping,SUMSTATS=${SUMSTATS},CUTOFF=${CUTOFF} \
--partition pfen1 --mem 40G ${CODEDIR}/slurm_polyfun_aggregate.sh
fi
done

