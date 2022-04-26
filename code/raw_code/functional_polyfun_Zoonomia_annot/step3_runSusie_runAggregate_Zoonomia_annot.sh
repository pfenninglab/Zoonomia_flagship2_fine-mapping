SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
CACHEDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/polyfun/LD_cache
DATADIR=${SETWD}/data/raw_data/functional_polyfun_Zoonomia_annot
CODEDIR=${SETWD}/code/raw_code/functional_polyfun_Zoonomia_annot
POLYFUNDIR='/home/bnphan/src/polyfun'
CUTOFF=5e-8; cd $CODEDIR;


#######################################
## split finemapping job by chromosome 
for ID in {1..34}; do
if [[ $ID -lt 27 ]]; then
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Loh_2018"
else
PREFIX=$(awk -F'\t' -v IND=${ID} 'FNR == IND + 1 {print $2}' ${SETWD}/data/tidy_data/tables/readme_ukbb_gwas.tsv)"-Gazal_2022"
fi
OUTDIR=${DATADIR}/${PREFIX}/susie
# 23G good for blocks w/ < 15k SNPs, 47G good for blocks <30k snps
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_Zoonomia_annot_aggregate.txt.gz && -f ${OUTDIR}/polyfun_all_jobs_22.txt ]]; then
echo "Submitting job for ${PREFIX} GWAS."
sbatch -p pool1,pfen1 --mem 45G --time 3-0:00:00 --array 1-22 \
--export=JobsFileName="${OUTDIR}/polyfun_all_jobs_@.txt" \
${SETWD}/code/raw_code/nonfunct_finemapping/slurm_finemap_byLine.sh
fi; done


###########################################
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
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_Zoonomia_annot_aggregate.txt.gz && -f ${OUTDIR}/polyfun_all_jobs_22.txt ]]; then
echo "Submitting job for ${PREFIX} GWAS."; 
sbatch -p pfen_bigmem,pool3-bigmem,pfen3 --mem 220G --time 4-4:00:00 --export=JobsFileName=${JobsFileName} \
${SETWD}/code/raw_code/nonfunct_finemapping/slurm_finemap_byLine.sh
fi; done


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
CUTOFF=5e-8
SUMSTATS=${SETWD}/data/tidy_data/polyfun/munged/${PREFIX}.parquet
OUTDIR=${DATADIR}/${PREFIX}/susie
if [[ ! -f ${DATADIR}/${PREFIX}/${PREFIX}_Zoonomia_annot_top_annot.txt.gz && -f ${OUTDIR}/polyfun_all_jobs_22.txt ]]; then
echo "Aggregating results from ${PREFIX} GWAS with P < ${CUTOFF} cutoff for loci."
sbatch --export=OUTDIR=${OUTDIR},PREFIX=${PREFIX}_Zoonomia_annot,SUMSTATS=${SUMSTATS},CUTOFF=${CUTOFF} \
--partition pfen1 --mem 23G ${SETWD}/code/raw_code/nonfunct_finemapping/slurm_polyfun_aggregate.sh
fi; 
done

