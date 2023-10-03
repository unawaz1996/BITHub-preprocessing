library(recount3)
library(magrittr)
library(tibble)
library(reshape2)
library(SummarizedExperiment)
library(corrplot)
library(dplyr)
library(ggvenn)
library(pander)
library(gridExtra)
library(variancePartition)
library(DT)
library(EnsDb.Hsapiens.v86)
library(singscore)
library(ggplot2)
library(UpSetR)
library(patchwork)
library(pheatmap)
library(ggpubr)
library(variancePartition)


brain.dir =file.path("/home/neuro/Documents/BrainData/Bulk")
e <- list()

e$Raw$GTEx <- read.csv(file.path(brain.dir, "/GTEx/Formatted/GTEx-exp.csv"),
                       check.names = FALSE, row.names = 1)
e$Raw$PE <- read.csv(file.path(brain.dir, "/PsychEncode/Formatted/PsychEncode-exp.csv"),
                     row.names =1, check.names = FALSE)
e$Raw$BSpan <- read.csv(file.path(brain.dir, "/BrainSpan/Formatted/BrainSpan-exp.csv"),
                        check.names = FALSE, row.names = 1)
e$Raw$BSeq <- read.csv(file.path(brain.dir,"/Brainseq/Formatted/BrainSeq-exp.csv"),
                       check.names = FALSE, row.names = 1)

e$Raw$HDBR <- read.csv(file.path(brain.dir, "/HDBR/Formatted/HDBR-exp.csv"),
                       check.names = FALSE, row.names = 1)  %>%
  column_to_rownames("EnsemblID")


md = list()
md$GTEx <- read.csv(file.path(brain.dir, "/GTEx/Formatted/GTEx-metadata.csv"),
                    check.names = FALSE, row.names = 1)
md$PE <- read.csv(file.path(brain.dir,"/PsychEncode/Formatted/PsychEncode-metadata.csv"),
                  check.names = FALSE, row.names=1)
md$BSeq <- read.csv(file.path(brain.dir,"/Brainseq/Formatted/BrainSeq-metadata.csv"),
                    check.names = FALSE)
md$HDBR <- read.csv(file.path(brain.dir, "/HDBR/Formatted/HDBR-metadata.csv"),
                    check.names = FALSE, row.names = 1)

md$BSpan <- read.csv(file.path(brain.dir,
                               "/BrainSpan/Formatted/BrainSpan-metadata.csv"),
                     check.names = FALSE)


### expression thresholding


gtex <- thresh(e$Raw$GTEx)
span <- thresh(e$Raw$BSpan)
seq <- thresh(e$Raw$BSeq)
pe <- thresh(e$Raw$PE)
hdbr <- thresh(e$Raw$HDBR)

## BrainSeq

varPar.bseq <- fitExtractVarPartModel(seq, form.bseq, md$BSeq)
varPar.bseq = sortCols(varPar.bseq)
plot_varPart = plotVarPart(varPar.bseq)

ggsave(plot_varPart,
       filename = here::here("output/Thesis_plots/Variance/PCA/BrainSeq/BrainSeq_varPart_full.svg"),
       height = 6,
       width = 8,
       units = "in")


write.csv(varPar.bseq, file = "/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted/BrainSeq-varPart.csv")

## BrainSpan

form.bspan <- ~ AgeNumeric  + (1|Sex) + (1|Period) + (1|Regions) + RIN +  pH + PMI + DonorID + Dissectionscore
varPar.bspan <- fitExtractVarPartModel(span,form.bspan ,md$BSpan)
varPar.bspan = sortCols(varPar.bspan)
plot_varPart = plotVarPart(varPar.bspan)

ggsave(plot_varPart,
       filename = here::here("output/Thesis_plots/Variance/PCA/BrainSpan/BrainSpan_varPart_full.svg"),
       height = 6,
       width = 8,
       units = "in")


write.csv(varPar.bspan, file="/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted/BrainSpan-varPart.csv")

## GTEx

gtex.form <- ~ TotalNReads + rRNA_ratio + (1|TypeofBatch)  + (1|AgeInterval) + (1|Sex) + (1|Regions) + RIN + HardyScale + PMI + MappingRate + PMI + ExonicRate

varPar.gtex <- fitExtractVarPartModel(gtex,gtex.form ,md$GTEx)
varPar.gtex = sortCols(varPar.gtex)
# varPar.gtex = read.csv("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-varPart.csv",
#                        header=TRUE,
#                        row.names = 1)
plot_varPart = plotVarPart(varPar.gtex)


ggsave(plot_varPart,
       filename = here::here("output/Thesis_plots/Variance/PCA/GTEx/GTEx_varPart_full.svg"),
       height = 6,
       width = 8,
       units = "in")

write.csv(varPar.gtex,"/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted/GTEx-varPart.csv")

##HDBR
form.hdbr <- ~ (1|AgeInterval) + (1|Sex) + (1|Hemisphere) + (1|SequencingBatch)+  TotalNReads + MappingRate + (1|Regions) + PMI + DonorID + IntronicRate
varPar.hdbr <- fitExtractVarPartModel(hdbr,form.hdbr ,md$HDBR)
varPar.hdbr = sortCols(varPar.hdbr)

plot_varPart = plotVarPart(varPar.hdbr)

ggsave(plot_varPart,
       filename = here::here("output/Thesis_plots/Variance/PCA/HDBR/HDBR_varPart_full.svg"),
       height = 6,
       width = 8,
       units = "in")

write.csv(varPar.hdbr , file = "/home/neuro/Documents/BrainData/Bulk/HDBR/Formatted/HDBR-varPart.csv")

## PsychEncode


md$PE= md$PE[!duplicated(md$PE[,c('SampleID')]),]
varPar.pe <- fitExtractVarPartModel(pe, pe.form, md$PE)

varPar.pe = sortCols(varPar.pe)
plot_varPart = plotVarPart(varPar.pe)
plot_varPart

ggsave(plot_varPart,
       filename = here::here("output/Thesis_plots/Variance/PCA/PsychEncode/PE_varPart_full.svg"),
       height = 6,
       width = 8,
       units = "in")

write.csv(varPar.pe,"/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted/PsychEncode-varPart.csv")
