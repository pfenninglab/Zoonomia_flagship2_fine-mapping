SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
DATADIR=${SETWD}/data/raw_data/functional_polyfun_baseline-LF2.2.UKB
CODEDIR=${SETWD}/code/raw_code/functional_polyfun_baseline-LF2.2.UKB
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
# 23G good for blocks w/ -lt 15k SNPs, 47G good for blocks -lt30k snps
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_baseline-LF2.2.UKB_aggregate.txt.gz ]]; then
echo "Submitting job for ${PREFIX} GWAS."; OUTDIR=${DATADIR}/${PREFIX}/susie
sbatch -p pool1,pfen1 --mem 45G --time 3-0:00:00 --array 1-22 \
	--export=JobsFileName="${OUTDIR}/polyfun_all_jobs_@.txt" \
	${SETWD}/code/raw_code/nonfunct_finemapping/slurm_finemap_byLine.sh
fi
done


##########################################
## combine all the jobs together per trait
## catch the jobs that need a lot of RAM
for ID in {1..34}; do
if [[ $ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
OUTDIR=${DATADIR}/${PREFIX}/susie
## merge all jobs together
JobsFileName=${OUTDIR}/polyfun_all_jobs.txt
if [[ ! -f $JobsFileName && -f ${OUTDIR}/polyfun_all_jobs_22.txt || $(wc -l $JobsFileName | cut -d ' ' -f1 ) == 0  ]]; then cat ${OUTDIR}/polyfun_all_jobs_*.txt > $JobsFileName; fi
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_baseline-LF2.2.UKB_aggregate.txt.gz ]]; then
echo "Submitting job for ${PREFIX} GWAS."
sbatch -p pool3-bigmem,pfen_bigmem,pfen1,pfen3 --mem 120G --time 4-4:00:00 --export=JobsFileName=${JobsFileName} \
	${SETWD}/code/raw_code/nonfunct_finemapping/slurm_finemap_byLine.sh
fi
done


##################################################
## submit job to aggregate fine-mapping when done
## catch the jobs that need a lot of RAM
for ID in {1..34}; do
if [[ $ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
N=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $4}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)
CUTOFF=5e-8; OUTDIR=${DATADIR}/${PREFIX}/susie
SUMSTATS=${SETWD}/data/tidy_data/polyfun/munged/${PREFIX}.parquet
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_baseline-LF2.2.UKB_top_annot.txt.gz ]]; then
echo "Aggregating results from ${PREFIX} GWAS with P -lt ${CUTOFF} cutoff for loci."
sbatch --export=OUTDIR=${OUTDIR},PREFIX=${PREFIX}_baseline-LF2.2.UKB,SUMSTATS=${SUMSTATS},CUTOFF=${CUTOFF} \
--partition pool1 --time 3-00:00:00 --mem 23G ${SETWD}/code/raw_code/nonfunct_finemapping/slurm_polyfun_aggregate.sh
fi
done

