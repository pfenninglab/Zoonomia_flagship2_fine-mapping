#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --job-name=run_annot
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --error=logs/run_annot_%A_%a.out.txt
#SBATCH --output=logs/run_annot_%A_%a.out.txt
#SBATCH --array=1-22

######################### USAGE #############################
#   sbatch /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh \
#    -i /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/foreground/mouse2human_tissue/mouse_atac-seq_peaks_upregulated_CTX_20190929_halLifted2hg38_orthologs.bed \
#    -n test -o /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/annotations/test
ulimit -u 8192;
function usage()
{
    echo "annotate_bed_LDSC_hg38.sh takes a bedfile in hg38 coordinates and annotates it with LDSC baseline files and HapMap3 SNPs for LDSC cell-type specific enrichments"
    echo "Uses slurm array jobs to submit the jobs to process for all 22 autosomal chromosomes"
    echo "Make sure `make_annot.py` and `ldsc.py` are in your ~/.bashrc path."
    echo "Example run call below:"
    echo ""
    echo "sbatch annotate_bed_LDSC_hg38.sh -i myBedFile_hg38.bed -n annotationName -o /dir/to/annotations/"
    echo "[-h]--help"
    echo "[-i]--input-bed-file      = BED-FILE"
    echo "[-w]--input-bigwig-file   = BIGWIG-FILE"
    echo "[-n]--name-base           = NAME-BASE"
    echo "[-o]--output-dir          = OUTDIR"
    echo ""
    echo "Note 1: install pyBigWig if in the ldsc conda envir if providing bigwig for continuous annotations:"
    echo "  conda activate ldsc"
    echo "  conda install pybigwig -c conda-forge -c bioconda --yes"
    echo ""
    echo "Note 2: use the '@' symbol where the numeric chromosome number is replaced if bigwig is split by chr."
    echo " e.g. -w myBigWig.@.bw  will stand in for {myBigWig.1.bw, myBigWig.2.bw,...}"
}

# read in command line arguments
while [ "$1" != "" ]; do
    case $1 in
        -i | --input-bed-file ) shift
                                BEDFILE=$1
                                ;;
        -n | --name-base )    	shift
                                NAME=$1
                                ;;
        -o | --output-dir )     shift
                                OUTDIR=$1
								;;
        -h | --help )           usage
                                exit 1
                                ;;
        *)
            usage
            exit 1
    esac
    shift

done

# make output directory if does not exist
source ~/.bashrc; source activate ldsc

# change the @ symbol to current chromosome number
mkdir -p ${OUTDIR}

# check to see if bedfile exists
if [ $(ls ${BEDFILE}| wc -l) == 0 ]
	then
    echo "Input bedfile ${BEDFILE} does not exist"
 	exit 1
elif [ $(ls ${BWFILE}| wc -l) == 0 ]; then
    echo "Input bigwig file ${BWFILE} does not exist."
    exit 1
fi

###############################################
### annotate the 1000G files with bed peaks and/or bigwig file
ANNOTFILE=${OUTDIR}/${NAME}.${SLURM_ARRAY_TASK_ID}.annot.gz
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
if [ "${BWFILE}" == "" ]; then
    echo "Partitioning baseline LD with ${NAME} bed regions"
    make_annot.py --bed-file ${BEDFILE} --annot-file ${ANNOTFILE} \
        --bimfile ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/plink_files/1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID}.bim
else
    echo "Partitioning with 1000G LD scores on HapMap3 SNPs with values from ${BWFILE} overlapping ${NAME} bed regions."
    ${GWASDIR}/scripts/make_annot_continuous.py --annot-name ${NAME} \
    --bigwig-file ${BWFILE} --bed-file ${BEDFILE} --annot-file ${ANNOTFILE} \
    --bimfile ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/plink_files/1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID}.bim
fi

##########################################
###  Compute LD scores for this annotation 
echo "Computing partitioned LD scores for ${NAME}."
ldsc.py --l2 --bfile ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/plink_files/1000G.EUR.hg38.${SLURM_ARRAY_TASK_ID} \
	--ld-wind-cm 1 --thin-annot --print-snps ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/listHM3.noMHC.txt \
	--annot ${ANNOTFILE} --out ${OUTDIR}/${NAME}.${SLURM_ARRAY_TASK_ID}

echo "Done."