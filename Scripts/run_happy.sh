#!/bin/bash

##########################################################################
#
# Platform: University of Sydney Artemis HPC
# Usage: qsub run_happy.pbs (or bash run_happy.sh)
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
#PBS -N 
#PBS -l select=1:ncpus=3:mem=10GB
#PBS -l walltime=03:00:00
#PBS -m e
#PBS -M 
#PBS -q defaultQ
#PBS -W umask=022
#PBS -e ./Logs/happy.e
#PBS -o ./Logs/happy.o

# load modules
module load hap/0.3.14
module load bcftools/1.14
module load htslib/1.14

# define variables
sample=
truth=
ref=
confidence=
out=${sample}.happy

# make logs directory in current directory
mkdir -p ./Logs

# remove INFO,FORMAT fields not used by hap.py from query vcf to avoid errors
echo "CLEANING UP QUERY VCF FOR HAPPY..."

bcftools annotate -Oz -x INFO,FORMAT ${samplevcfdir}/${sample}.vcf.gz -o - | \
        bcftools sort -Oz - -o ${samplevcfdir}/${sample}_geno_sorted.vcf.gz

# index sorted query vcf
echo "INDEXING QUERY VCF..."

tabix -f ${samplevcfdir}/${sample}_geno_sorted.vcf.gz

# remove overlapping variants in truth vcf - these throw an error and hap.py fails
echo "PREPARING TRUTH VCF TO AVOID ERRORS..."

bcftools view -e "%FILTER='SiteConflict'" -Oz ${work}/${truth}.vcf.gz -o ${work}/${truth}_clean.vcf.gz

# index cleaned truth vcf
tabix ${work}/${truth}_clean.vcf.gz

echo "RUNNING HAP.PY VCF COMPARISON TOOL..."

# run hap.py using cleaned truth vcf and sample vcf
hap.py ${work}/${truth}_clean.vcf.gz \
        ${samplevcfdir}/${sample}_geno_sorted.vcf.gz \
        --report-prefix ${out} \
        -r ${ref} \
        -L --preprocess-truth --usefiltered-truth --fixchr \
        --roc QUAL \
        --engine xcmp \
        --verbose
