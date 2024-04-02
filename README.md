# BITHub
Brain Integrative Transcriptome Hub (BITHub) is a web resource that aggregates gene-expression data from the human brain across multiple consortia, and allows direct comparison of gene expression in an interactive manner. This repository contains the code used for preprocessing transcriptomic data used for BITHub. 

<br>

Table of Contents:

- [Datasets](#datasets)
  * [Data collection](#data-collection)
- [Preprocess](#preprocessing)
  * [Metadata annotation](#metadata-annotation)
  * [Determining drivers of variation](#determining-drivers-of-variation)
  * [Normalization](#normalization)


# Datasets 


## Data collection 

Both processed bulk and single-nucleus RNA-seq human brain transcriptomic datasets were retrieved from their respective portals as highlighted in Table 1. 

| Dataset   | Description | nSamples | Original file |
| ------------- | ------------- | ------------- | ------------- |
| BrainSeq  | RNA-seq data of the human postmortem brain including hippocampus and dorsolateral prefrontal cortex. Collado-Torres et al used RiboZero libraries on 900 tissue samples from 551 individuals (including 286 with schizophrenia). Prenatal  (age < 0; range, 14 to 22 post-conception weeks) and postnatal (age ≥ 18 years; range, 18 to 96 years) samples were used in this work. Resource: [BrainSeq Phase II](https://eqtl.brainseq.org/phase2/) <br><br>| 900 |  [Expression matrix and metadata](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata) |
| BrainSpan  | Samples collected and analyzed by Kang et al across multiple brain structures including 11 neocortical areas, cerebellar cortex, mediodorsal nucleus of the thalamus, striatum, amygdala, and hippocampus. Samples included prenatal (age < 0; range, 8 to 38 post-conception weeks) and postnatal (age ≥ 4 mos ; range, 4 mos to 41 years) phenotypes of the normal human brain   | 524  | [BrainSpan Developmental Atlas](https://www.brainspan.org/static/download.html) <br><br> [Expression matrix and metadata](http://www.brainspan.org/api/v2/well_known_file_download/267666525) <br><br> Additional metadata information <br> [Allen Brain Atlas](https://help.brain-map.org/download/attachments/3506181/Human_Brain_Seq_Stages.pdf?)  <br><br> mRIN: <br> Feng et al (2015)|
| GTEx  | The Genotype-Tissue Expression database contains 2,642 samples of the human postmortem brain in postnatal ages (age <20; range 20 to 79 years across 13 brain regions. All samples have been collected from non-diseased individuals   | 2642  | [GTEx v8](https://gtexportal.org/home/datasets) <br><br> [Gene TPMs](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz) <br> <br> Metadata files <br> [Phenotype Attributes](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt) <br> [Sample Attributes](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt) |
| PsychEncode | The PsychEncode dataset contains data of the dorsolateral prefrontal cortex from human postmortem tissues from prenatal (age <0; range 4 to 40 pcw) and postnatal samples. Samples include controls and individuals with diagnosis of Bipolar Disorder, ASD, Schizophrenia and Affective disorder.   | 1866  | [PsychEncode Resource](http://resource.psychencode.org) <br><br> [Expression matrix](http://resource.psychencode.org/Datasets/Derived/DER-02_PEC_Gene_expression_matrix_TPM.txt)  <br><br> Metadata file  <br> *Access from Synapse* |
| Human Cell Atlas  | Content Cell  | 32,749 | Content Cell  |
| Velmeshev et al  | Velmeshev et al generated single-nuclei from 48 post-mortem tissue samples from the prefrontal cortex, anterior cingulate and insular cortical regions. Donors included 16 control subjects and 11 patients with ASD. All samples are postnatal  | 81,216 | [Cells UCSC](https://cells.ucsc.edu/?ds=autism)  <br> <br> Matrix: <br> [exprMatrix.tsv.gz](https://cells.ucsc.edu/autism/exprMatrix.tsv.gz) <br> Values in matrix are: 10x UMI counts from cellranger, log2-transformed <br><br> Raw count matrix: <br> [rawMatrix.zip](https://cells.ucsc.edu/autism/rawMatrix.zip) |

# Preprocessing 

All scripts for pre-processing data are loading in this Github repository. In order to use the code for preprocessing, please check the workflowR document: 
The `analysis/` subfolder within this project contains detailed code for each step of the preprocessing and analysis for files for BITHub, including code for the analysis for the thesis chapter. Prior to rerunning the Rmarkdown files, please change the location of input and output directories in the 'code/directories.R' file. 


The following Rmarkdown files are provided: 

| File name  | Description | output file | 
| ------------- | ------------- | ------------- | 
`01-bulk-preprocess-data.Rmd` | This file contains the detailed analysis of preprocessing metadata and expression matrices of the bulk RNA-seq datasets listed in Table 1. The input of each file is the expression matrix as provided by each given dataset / database, and original metadata. The output is a harmonised metadata, whereby consistent naming conventions for metadata variables across datasets are adopted. In addition, consistent brain related nonemclature and developmental stages with respect to sample information are also provided. Full details are [Metadata annotation](#metadata-annotation) | |
|`02-snRNA-seq-preprocess.Rmd` |||
|`03-metadata-attributes.Rmd` |||
||||
||||


## Metadata annotation

As the metadata annotation was heterogeneous across the datasets, rigorous harmonization was performed. For each dataset, columns specifying `Age Intervals`, `Regions`, `Diagnosis` and `Period` were also added.

<br> 

*Developmental Ages*
<br>
Samples were binned into age intervals that were used to define developmental stages. For all samples < 20 years old, the binning was performed based on the BrainSpan Technical White Paper (Kang el al, 2011), whereas samples  20 years were binned in 10 year intervals. 

To allow comparison on a consistent scale, all ages were converted to years (numeric age). 
For prenatal ages (labeled -pcw): 

```
Numeric age =  -(40 - pcw) 52 
```

where pcw is the age in post-conception weeks, 40 is the total number of prenatal weeks, and 52 denotes the total number of weeks in a year. 

For ages labelled in months (labelled mos): 

```
Numeric age = mos / 12
```

where 12 represents the total number of months in a year. 

Prenatal or postnatal tags were assigned to samples depending on their numeric age where numeric age < 0 was labeled as prenatal and numeric age < 0 as postnatal.


Ontology and nomenclature of brain regions <br>

The brain structures were divided into 4 main categories (regions):` Cortex`, `Subcortex`, `Cerebellum` and `Spinal Cord`. 


## Determining drivers of variation

`variancePartition` was used for mixed linear analysis to estimate the proportion of variance explained by the selected covariates on each gene.  Highly correlated covariates cannot be included in the model, and so as a result covariates that were not strongly correlated for each dataset were selected to the the varianceParititon analysis on filtered genes. For filtering, the expression cut-off was selected at 1 RPKM/TPM/CPM  in at least 10% of the samples

## Normalization 
To allow direct comparison of datasets with different normalizations, datasets have `z-score transformed mean log2 expression values`. 
