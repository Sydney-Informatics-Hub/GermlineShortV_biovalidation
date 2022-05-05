# GermlineShortV_biovalidation

 - [Description](#description)
  - [Diagram](#diagram)
  - [User guide](#user-guide)
      - [Quick start guide](#quick-start-guide)
  - [Benchmarking](#benchmarking)
  - [Workflow summaries](#workflow-summaries)
      - [Metadata](#metadata)
      - [Component tools](#component-tools)
      - [Required (minimum)
        inputs/parameters](#required-minimum-inputsparameters)  
        [Preparing your own input files](#preparing-input-files)
  - [Additional notes](#additional-notes)
      - [Understanding your outputs](#understanding-your-outputs)  
      - [Performance metrics explained](#performance-metrics-explained)   
  - [Help/FAQ/Troubleshooting](#helpfaqtroubleshooting)
  - [Acknowledgements/citations/credits](#acknowledgementscitationscredits)

## Description 
Population-scale WGS cohorts are essential resources for genetic analyses including heritable diseases, evolutionary genomics, conservation biology, and population genomics. Processing raw reads into analysis-ready variants remains challenging. Various mapping and variant calling pipelines have been made publicly available in recent decades. Designing a mapping and variant calling pipeline to meet your needs is dependent on the compute infrastructure you’re working on, the types of variants you’re primarily interested in, and the sequencing technology you use to generate raw sequencing data. Keep in mind that the tools you use to build your pipeline can affect variant calling accuracy. Further, optimisation and customisation of these tools’ commands can also affect their performance. Best-practice recommendations for variant calling pipelines vary dramatically between species and research questions, depending on the availability of genomic resources for the population of interest, genome structure, and clinical relevance of the resulting variant dataset. It is important to not only design a robust variant calling pipeline but also fine-tune it to achieve optimal performance for your dataset and research question. 

There are various measurements that you can apply to evaluate the biological accuracy of your germline variant calling pipeline. Currently, no best practice methods for interrogating joint-called variant sets exist in the literature. A number of publicly available, human ‘gold standard’ truth datasets including Platinum Genomes and Genome in a Bottle (GIAB) are useful for benchmarking across high confidence regions of the genome and evaluating the recall and precision of the pipeline. We recommend individuals working with human datasets benchmark their germline variant calling pipelines using one of these datasets. Unfortunately, these resources are not typically available for non-human organisms. 

Here, we present protocols for benchmarking and validating germline short variant (SNVs and indels) datasets using a combination of methods that can capture the quality of your variant sets for human, non-human model, and non-model organisms. The process you can apply will depend on the organism you’re working with and the genomic resources available to that organism. 

## Diagram 

<p align="center"> 
<img src="https://github.com/Sydney-Informatics-Hub/GermlineShortV_biovalidation/blob/main/Benchmarking%20and%20validation%20protocol.png" width="70%" height="70%">  
</p> 

## User guide 
###  Quick start guide 

These bash scripts were written for the University of Sydney’s high performance computer, Artemis. They can be run on the command line or submitted as PBS jobs. These scripts assume your input is a gzipped multi-sample (cohort) VCF file. Before running, edit the PBS project directive and define the variables at the top of the script. All software used in this protocol is installed on Artemis- to use alternate versions or run on a different compute infrastructure, edit the modules according to your needs.  

#### Human datasets 
For human datasets, we recommend you benchmark your germline variant calling pipeline using a gold standard dataset such as Platinum Genomes. Raw sequence data in FASTQ format for these datasets can be downloaded along with their high confidence variant calls and regions from public repositories. See [Preparing input files]() for more information on how to download and prepare these files.    

##### 1. Collect vcf summary metrics  
Edit the PBS -P directive and variables for your dataset in `vcfstat.sh`. Then run script with: 

```
qsub vcfstat.sh (or bash vcfstat.sh)
```
This will produce summary and quality metrics reports and plots for your cohort. It will also produce summary and detail files for known variant representation. BCFtools stats plots will be housed in a directory labelled `${cohort}_vcfplots`. 

##### 2. Biological benchmarking using a truth set  

Edit the PBS -P directive and variables for your files. Then run script with:  

```
qsub run_happy.sh
```
This script will subset your multi-sample VCF into individual samples, prepare them for hap.py, and output a number of files including summary metrics (including recall, precision and F1-score) and ROC count files that can be used to produce ROC curves, separately for SNVs and indels. See the [hap.py user guide](https://github.com/Illumina/hap.py/blob/master/doc/happy.md) for more information on how to interpret hap.py output. ROC curves of Hap.py runs can be plotted using the script [rocplot.Rscript](https://github.com/Illumina/hap.py/blob/master/src/R/rocplot.Rscript).   

#### Non-human model organism datasets

##### 1. Collect vcf summary metrics  
Edit the PBS -P directive and variables for your dataset in `vcfstat.sh`. We recommend you use the set of known variants used for base quality score recalibration to validate population level variants. If you used trio data, unhash the Mendelian error command within the script. Then run script with: 

```
qsub vcfstat.sh (or bash vcfstat.sh)
```
This will produce summary and quality metrics reports and plots for your cohort. It will also produce summary and detail files for known variant representation. BCFtools stats plots will be housed in a directory labelled `${cohort}_vcfplots`.  
#### Non-model organism datasets 

##### 1. Collect vcf summary metrics  

Edit the PBS -P directive and variables for your dataset in `vcfstat_nonmodel.sh`. Then run script with: 

```
qsub vcfstat_nonmodel.sh (or bash vcfstat_nonmodel.sh)
```

This will produce summary and quality metrics reports and plots for your cohort. It will also produce summary and detail files for known variant representation. BCFtools stats plots will be housed in a directory labelled `${cohort}_vcfplots`. 

## Benchmarking 
Coming soon!  

## Workflow summaries 
### Metadata 
|metadata field     | workflow_name / workflow_version  |
|-------------------|:---------------------------------:|
|Version            | 1.0                 |
|Maturity           | stable                            |
|Creators           | Georgie Samaha, Tracy Chew, Cali Willet                 |
|Source             | NA                                |
|License            | NA                                |
|Workflow manager   | NA                          |
|Container          | None                              |
|Install method     | Manual                            |
|GitHub             | NA                                |
|bio.tools 	        | NA                                |
|BioContainers      | NA                                | 
|bioconda           | NA                                |

### Component tools 

bcftools/1.14  
htslib/1.14  
python/3.8.2  
R/4.1.1  
hap.py/0.3.14  

### Required (minimum) inputs/parameters 

- Multi-sample or single sample VCF file (VCF.gz format)
- List of sample IDs that match the VCF (.txt format)
- Known variant dataset (VCF format. Human and non-human model organisms only)
- Pedigree file (format: mother,father,offspring. Trios or Platinum Genomes only)
- Truth set variant calls (VCF.gz format. Human, Platinum Genomes only)
- High confidence call regions (BED format. Human, Platinum Genomes only)

### Preparing input files 

#### Gold standard variant truth sets  

The benchmarking protocol for human datasets assumes you have performed mapping and germline variant calling on a gold standard truth set. These datasets contain millions of variants that have been confirmed using orthologous technologies [Eberle et al. 2017](https://doi.org/10.1101/gr.210500.116).   

We recommend you use the Platinum Genomes dataset for benchmarking germline variant calling pipelines that include joint genotyping of multiple samples. Six members, comprising two trios, of the Platinum Genomes dataset can be downloaded from the Illumina BaseSpace Sequence Hub, the ENA, or dbGaP. The Platinum Genomes dataset contains multiple files including the following files you will need for running `run_happy.sh`: 
- Paired-end FASTQ files for each sample
- High-confidence germline variant VCF files for each sample
- High-confidence genomic regions (BED format)

Currently, these files are available for Hg19 (GRCh37) and Hg38 (GRCh38) . Links to raw data are [here](https://github.com/Illumina/PlatinumGenomes). BaseSpace offers a command line tool for downloading files, see [here](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-examples) for instructions. 

#### Providing your own ‘truth set’ 
*A word of caution*- testing the performance of your pipeline using a truth set is only intended to estimate the overall quality of your pipeline and detect any potential sources of error in your method. It is not intended to test the truthfulness of your variant set. See [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531572-Evaluating-the-quality-of-a-germline-short-variant-callset) for further discussion of the assumptions we make about truth sets. Most non-human organisms do not have access to gold standard truth set resources like the Platinum Genomes dataset. However there are a few alternative options you could try: 
 - Genotyping arrays: if you have genotyping data for the same samples you tested your germline variant calling pipeline with, you can reformat these to VCF using a tool like [PLINK’s recode](https://www.cog-genomics.org/plink/1.9/data#recode) and use it as a truth set. 
 - Known variant datasets: if your organism of interest has a set of known population-level variants you can use these as a truth-set. Just remember that these variants might not always be validated (i.e. dbSNP). 

Using this method you will need to also provide your own high-confidence regions file in BED format. The location and size of these regions will depend on your dataset, organism, reference assembly and sequencing method. Typically these regions would exclude centromeres, telomeres and repetitive parts of the genome that are likely to complicate variant calling.   


## Additional notes 

Test data for Hap.py can be found [here](https://github.com/Illumina/hap.py/blob/master/doc/microbench.md)  

Instructions on how to install Hap.py can be found [here](https://github.com/Illumina/hap.py#installation)   

This warning may be thrown by Hap.py and can be ignored: `WARNING  No reference file found at default locations. You can set the environment variable 'HGREF' or 'HG19' to point to a suitable Fasta file.`  


### Understanding your outputs 
The following files will be produced and stored in your designated working directory. They will all be labelled with your specified cohort name.  

#### Variant based metrics 
Produced by BCFtools stats command. Output file:
- ${cohort}.bcftools.metrics  
- ${cohort}_bcftools.metrics_vcfstatplots (directory and files)  

#### Sample based metrics   
Produced by BCFtools smplstats and mendelian commands. Output files:
- ${cohort}.smplstats
- ${cohort}.smplstats.pdf
- ${cohort}.Mendelianerr

#### Known variant concordance 
Produced by GATK CollectVariantCallingMetrics command. Output files:
- ${cohort}.known.variant_calling_summary_metrics
- ${cohort}.known.variant_calling_detail_metrics

#### Biological validation using a truth set 
Produced by Hap.py. Output files:
- ${sample}.happy.metrics.json.gz
- ${sample}.happy.roc.all.csv.gz
- ${sample}.happy.roc.Locations.INDEL.csv.gz
- ${sample}.happy.roc.Locations.INDEL.PASS.csv.gz
- ${sample}.happy.roc.Locations.SNP.csv.gz
- ${sample}.happy.roc.Locations.SNP.PASS.csv.gz
- ${sample}.happy.roc.tsv
- ${sample}.happy.runinfo.json
- ${sample}.happy.summary.csv

### Performance metrics explained  

|Metric                                |Expected/ideal value                                |Tool           |Relevance                                                                                                      |
|--------------------------------------|----------------------------------------------------|---------------|---------------------------------------------------------------------------------------------------------------|
|Number of SNVs and indels (per sample)|Human WGS: ~4.4M, Human WES: ~41k, Species dependent|bcftools stats |Population, sequencing approach, and genomic region dependent. Alone, this metric cannot indicate data quality.|
|Indel length distribution             |Indel length range is 1-10,000bp.                   |bcftools stats |Increased length is conflated with reduced mapping quality. Distribution is dataset dependent. Recommend filtering for high quality.|
|Depth of coverage                     |Depends on the sequencing coverage of samples.      |bcftools stats |Dramatic deviation from expected distribution can indicate artifactual bias.                                   |
|Substitution type counts              |See TiTv ratio.                                     |bcftools stats |Twice as many possible transversions as transitions. See [here](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtu668)  |
|TiTv ratio (genome wide)              |For mammals: WGS: 2.0-2.1, WES: 3.0-3.3             |bcftools stats |Dramatic deviation from expected ratio can indicate artifactual bias. Typically elevated in coding regions where transversions are more likely to occur. |
|Base quality distribution             |Dataset dependent.                                  |bcftools stats |This will reflect the quality based filtering you performed. Dramatic deviation from expected ratio can indicate artifactual bias.|
|Indel ratio                           |Common: ~1.0, Rare: 0.2-0.5                         |GATK CollectVariantCallingMetrics|This should be evaluated after custom filtering variants for your needs. Dramatic deviation from expected ratio can indicate artifactual bias.|
|Het/hom(non-ref)                      |~2.0 assuming Hardy-Weinberg equilibrium.           |GATK CollectVariantCallingMetrics|Ancestry dependent, can vary dramatically. See [Wang et al. 2015](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtu668)|
|Mendelian error                       |0                                                   |BCFtools +mendelian|Mendelian inheritance errors are likely erroneous genotype calls. See [Pilipenko et al. 2014](https://dx.doi.org/10.1186%2F1753-6561-8-S1-S21)|
|True positives                        |Dataset dependent.                                  |Hap.py         |Number of query variants that are present in the truth set.                                                    |
|False negatives                       |Dataset dependent.                                  |Hap.py         |Number of variants in truth set, not present in query VCF.                                                     |
|False positives                       |Dataset dependent.                                  |Hap.py         |Number of variants in query VCF, not present in truth set.                                                     |
|Recall                                |1                                                   |Hap.py         |Absence of false negatives. See [Krusche et al. 2019](https://doi.org/10.1038/s41587-019-0054-x)               |
|Precision                             |1                                                   |Hap.py         |Absence of false positives. See [Krusche et al. 2019](https://doi.org/10.1038/s41587-019-0054-x)               |
|F1-score                              |1                                                   |Hap.py         |Harmonic mean of recall and precision. See [Krusche et al. 2019](https://doi.org/10.1038/s41587-019-0054-x)    |
|Genotype errors (FP.GT)               |Dataset dependent.                                  |Hap.py         |Number of query variants with incorrect genotype                                                               |

### Resources and references 

Eberle, M. A., Fritzilas, E., Krusche, P., Källberg, M., Moore, B. L., Bekritsky, M. A., Iqbal, Z., Chuang, H. Y., Humphray, S. J., Halpern, A. L., Kruglyak, S., Margulies, E. H., McVean, G., & Bentley, D. R. (2017). A reference data set of 5.4 million phased human variants validated by genetic inheritance from sequencing a three-generation 17-member pedigree. Genome research, 27(1), 157–164. https://doi.org/10.1101/gr.210500.116   

Koboldt, D.C. Best practises for variant calling in clinical sequencing. Genome Med 12, 91 (2020). https://doi.org/10.1186/s13073-020-00791-w  

Krusche, P., Trigg, L., Boutros, P.C. et al. Best practices for benchmarking germline small-variant calls in human genomes. Nat Biotechnol 37, 555–560 (2019). https://doi.org/10.1038/s41587-019-0054-x  

Marshall, C.R., Chowdhury, S., Taft, R.J. et al. Best practices for the analytical validation of clinical whole-genome sequencing intended for the diagnosis of germline disease. npj Genom. Med. 5, 47 (2020). https://doi.org/10.1038/s41525-020-00154-9   

Pilipenko, V.V., He, H., Kurowski, B.G. et al. Using Mendelian inheritance errors as quality control criteria in whole genome sequencing data set. BMC Proc 8, S21 (2014). https://doi.org/10.1186/1753-6561-8-S1-S21   

Wang, J., Raskin, J., Samuels, D., Shyr, Y., Guo, Y., Genome measures used for quality control are dependent on gene function and ancestry, Bioinformatics 31, 318–323 (2015)  https://doi.org/10.1093/bioinformatics/btu668  


## Help/FAQ/Troubleshooting

If Hap.py throws an error, search the [issues at Hap.py GitHub repository](https://github.com/Illumina/hap.py/issues) and attempt to resolve it before submitting an issue here.    

## Acknowledgements/citations/credits  

### Authors 
- Georgie Samaha (Sydney Informatics Hub)   
- Tracy Chew (Sydney Informatics Hub)  
- Cali Willet (Sydney Informatics Hub)  
- Nandan Deshpande (Sydney Informatics Hub)

Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:  

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia.  
