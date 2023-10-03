library(CePa)
library(EnsDb.Hsapiens.v86)
library(recount3)
library(DT)
library(dplyr)
library(tidyr)
source("functions.R")
library(pander)
library(gridExtra)
library(variancePartition)
library(corrplot)
library(edgeR)

### Load the data from GTEx 


gtex_trans = read.gct("/home/neuro/Documents/BrainData/Bulk/GTEx/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")
colnames(gtex_trans ) <- gsub("\\.", "-", colnames(gtex_trans))


dir = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/")
annot= read.csv("/home/neuro/Documents/Brain_integrative_transcriptome/BITHub-preprocessing/data/annotations/GTEx-metadata-annot.csv")
attributes = list.files(dir, full.names = TRUE, pattern = "\\SampleAttributesDS.txt") # Sample attributes contains sample level information
phenotype = list.files(dir, full.names = TRUE, pattern = "\\SubjectPhenotypesDS.txt")
md = read_tsv(attributes, col_types = c('.default' = 'c')) %>% 
    #filter(SMTS == 'Brain') %>% 
    mutate(SUBJID = sapply(str_split(SAMPID, pattern = "-"), function(x) paste(x[1:2], collapse = '-'))) %>%
    left_join(read_tsv(phenotype, col_types = c('.default' = 'c')))  %>% as.data.frame()
colnames(md) = annot$BITColumnName[match(colnames(md), annot$OriginalMetadataColumnName)]


gtex_trans = gtex_trans %>%
  as.data.frame() %>%
  rownames_to_column("tx_id") %>%
  mutate_at(.vars = "tx_id", .funs = gsub, pattern = "\\.[0-9]*$", replacement = "") %>% 
  column_to_rownames("tx_id") %>% 
  dplyr::select(contains(md$SampleID))

message("averaging GTEx per region now")
a = with(md, setNames(StructureAcronym, SampleID))[names(gtex_trans)[col(gtex_trans)]]
#a <- with(md, setNames(StructureAcronym, SampleID))[names(gtex_trans_update)[col(gtex_trans_update)]]
#average = tapply(unlist(gtex_trans), list(row(gtex_trans), a), mean)



message("saving file now")

save(a, file = "/home/neuro/Documents/Brain_integrative_transcriptome/BITHub-preprocessing/output/Results/transcript/transcript-gtex-average.Rda")