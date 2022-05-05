##########################################################################
#
# Platform: University of Sydney Artemis HPC
# Usage: this script is run by vcfstat.sh and vcfstat_nonmodel.sh
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

# this script is run by vcfstat.sh scripts.
# it plots summary info from $cohort.smplstats

## load libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)
library(ggsci)

args <- commandArgs(TRUE)

## Load smplstat file
samplelines <- read.delim(args[1], header=FALSE, comment.char="#", skip=29)

# collect only lines reporting sample-level info (FLT0)
FLT0<-samplelines[grep("FLT0", samplelines$V1),]

# split and define columns for sample info
colnames(FLT0)<-c("FLT", "IID", "TOTAL", "NON-REF", "HOM-REF", "HOM-ALT", "HET", "HEMI", "SNPs", "INDELs", "SINGLETONs", "MISS", "TRANSVERSIONS", "TRANSITIONS", "TITV")

## Make plots
# plot all samples total: snp vs indel counts for each sample
snpindelmelt<-melt(FLT0, id.vars= "IID", value.name="count", variable.name = "type", measure.vars = c("SNPs", "INDELs"))

snpsindels<-ggplot(snpindelmelt, aes(x=IID, y=count, fill=type))+
 geom_col(position='dodge')+
 theme_bw()+
 theme(panel.grid = element_blank(),
 axis.text.x = element_text(angle=90))+
 labs(title = "", x="", y="Variant count")+
 scale_fill_ucscgb()

# plot missing per sample
missing<-ggplot(FLT0, aes(x=IID,y=MISS, fill=IID))+
 geom_bar(stat = "identity")+
 theme_bw()+theme(panel.grid = element_blank(),
 axis.text.x = element_text(angle=90))+
 labs(title = "", x="", y="Number of missing sites")+
 scale_fill_ucscgb()

# plot transition to transversion counts for each sample
transmelt<-melt(FLT0, id.vars= "IID", value.name="count", variable.name = "type", measure.vars = c("TRANSVERSIONS", "TRANSITIONS"))

transcounts<-ggplot(transmelt, aes(x=IID, y=count, fill=type))+
 geom_col(position='dodge')+
 theme_bw()+theme(panel.grid = element_blank(),
 axis.text.x = element_text(angle=90))+
 labs(title = "", x="", y="Base substitution count")+
 scale_fill_ucscgb()

# plot all samples TITV
titv<-ggplot(FLT0, aes(x=IID,y=TITV, fill=IID))+
 geom_bar(stat = "identity")+
 geom_text(aes(label=TITV), vjust=0)+
 coord_cartesian(ylim = c(1.0, 2.0)) +
 theme_bw()+theme(panel.grid = element_blank(),
 axis.text.x = element_text(angle=90))+
 labs(title = "", x="", y="TiTv ratio")+
 scale_fill_ucscgb()

# arrange plots and save to PDF
smplstats <- ggarrange(snpsindels, missing, transcounts, titv, ncol = 2, nrow=2)
ggsave(file=paste0(args[1], ".pdf"), smplstats, width = 8, height = 8, dpi = 400, units = "in", device='pdf')
