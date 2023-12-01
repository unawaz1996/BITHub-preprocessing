library(CePa)
library(EnsDb.Hsapiens.v86)
library(recount3)
library(DT)
library(dplyr)
library(tidyr)
source("/home/neuro/Documents/Brain_integrative_transcriptome/BITHub-preprocessing/code/preprocess/functions.R")
library(pander)
library(gridExtra)
library(variancePartition)
library(corrplot)
library(edgeR)

### Load the data from GTEx


txdf = transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene = as.data.frame(txdf[,c("tx_id","gene_id", "tx_biotype")])

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


gtex_trans %<>%
  as.data.frame() %>%
rownames_to_column("tx_id") %>%
  mutate_at(.vars = "tx_id", .funs = gsub, pattern = "\\.[0-9]*$", replacement = "") %>% left_join(tx2gene, by="tx_id")

genes_in_gtex = unique(gtex_trans$gene_id)


tx_avg = list()
for (g in genes_in_gtex) {
  temp = gtex_trans %>%
    dplyr::filter(gene_id == g) %>%
    column_to_rownames("tx_id") %>%
    dplyr::select(-c(gene_id,tx_biotype)) %>%
    dplyr::select(contains(md$SampleID))
    message(paste("averaging GTEx per region now for", g))
    average = with(md, setNames(StructureAcronym, SampleID))[names(temp)[col(temp)]]

    average_exp = tapply(unlist(temp), list(row(temp), average ), mean)

    rownames(average_exp) = rownames(temp)
    tx_avg[[paste0(g)]] = average_exp

}


tx_avg_all = tx_avg %>%
  do.call(rbind, tx_avg)

write.csv(tx_avg_all, file = "/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-tx.csv")
#gtex_trans = gtex_trans %>%
#  as.data.frame() %>%
#  rownames_to_column("tx_id") %>%
#  mutate_at(.vars = "tx_id", .funs = gsub, pattern = "\\.[0-9]*$", replacement = "") %>%
#  column_to_rownames("tx_id") %>%
#  dplyr::select(contains(md$SampleID))


#a = with(md, setNames(StructureAcronym, SampleID))[names(gtex_trans)[col(gtex_trans)]]
#a <- with(md, setNames(StructureAcronym, SampleID))[names(gtex_trans_update)[col(gtex_trans_update)]]
#average = tapply(unlist(gtex_trans), list(row(gtex_trans), a), mean)



#message("saving file now")

#save(a, file = "/home/neuro/Documents/Brain_integrative_transcriptome/BITHub-preprocessing/output/Results/transcript/transcript-gtex-average.Rda")
