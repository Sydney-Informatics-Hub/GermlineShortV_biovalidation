#!/bin/bash

##########################################################################
#
# Platform: University of Sydney Artemis HPC
# Usage: qsub vcfstat.sh (or bash vcfstat.sh)
# Version: 1.0
#
# For more details see:
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the support provided by the Sydney Informatics Hub,
# a Core Research Facility of the University of Sydney. This research/project
# was undertaken with the assistance of resources and services from
# the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia funding.
#
##########################################################################

#PBS -P
#PBS -N vcf_summary
#PBS -l select=1:ncpus=1:mem=8GB
#PBS -l walltime=03:00:00
#PBS -m e
#PBS -M
#PBS -q defaultQ
#PBS -W umask=022
#PBS -e ./Logs/vcf_stats.e
#PBS -o ./Logs/vcf_stats.o

set -e

# define variables
vcfdir=/project/SIHsandbox/Platinumgenomes/VCFs #where are your vcfs stored
summariesdir=/project/SIHsandbox/Platinumgenomes/vcfsummaries #where do you want the summaries to be stored
pedigree= #whats the path to and name of your pedigree file
cohort=platinumgenomes_pawsey #whats the prefix of your cohort vcf
known= #what is the known variant database file you are using
mendelerr=${summariesdir}/${cohort}.Mendelianerr
bcftoolsmetrics=${summariesdir}/${cohort}.bcftools.metrics
popmetrics=${summariesdir}/${cohort}.known.metrics
smplstats=${summariesdir}/${cohort}.smplstats

# modules
module load bcftools/1.14
module load gatk/4.2.1.0
module load htslib/1.14
module load R/4.1.1

# make a directory to store PBS logs
mkdir -p ./Logs
mkdir -p $summariesdir

# cohort vcf QC and summary metrics - will also produce summary plots
bcftools stats ${vcfdir}/${cohort}.vcf.gz \
        > ${bcftoolsmetrics}
plot-vcfstats -p ${bcftoolsmetrics}_vcfstatsplots \
        ${bcftoolsmetrics}

# sample level vcf QC and summary metrics
bcftools plugin smpl-stats \
                --output ${smplstats} \
                ${vcfdir}/${cohort}.vcf.gz

Rscript smplstats.R ${smplstats}

# compare vcf with known population callset
gatk CollectVariantCallingMetrics \
        -I ${vcfdir}/${cohort}.vcf.gz \
        --DBSNP ${known} \
        -O ${popmetrics}

# Mendelian errors, unhash for trios vcf only
bcftools +mendelian ${vcfdir}/${cohort}.vcf.gz \
        -T ${pedigree} \
        -m c  > ${mendelerr}

# print summary to stdout
if [[ -d ${bcftoolsmetrics}_vcfstatsplots && -f ${smplstats} && -f ${popmetrics}.variant_calling_summary_metrics && -f ${popmetrics}.variant_calling_detail_metrics ]]; then
        printf "\nVCF SUMMARY STATS FINISHED\n\n"

printf "see ${bcftoolsmetrics}_vcfstatsplots directory for summary plots\n\n"
printf "see ${smplstats} for individual sample stats\n\n"
printf "see ${popmetrics} for known variant summary stats\n\n"

else test -d ${bcftoolsmetrics}_vcfstatsplots && echo ${bcftoolsmetrics}_vcfstatsplots was successfully produced || echo something went wrong, ${bcftoolsmetrics}_vcfstatsplots was not created
     test -f ${smplstats} && echo ${smplstats} was successfully produced || echo something went wrong, ${smplstats} was not created
     test -f ${popmetrics}.variant_calling_summary_metrics && echo ${popmetrics}.variant_calling_summary_metrics was successfully produced || echo something went wrong, ${popmetrics}.variant_calling_summary_metrics was not created
     test -f ${popmetrics}.variant_calling_summary_metrics && echo ${popmetrics}.variant_calling_detail_metrics was successfully produced || echo something went wrong, ${popmetrics}.variant_calling_detail_metrics was not created

fi
