# List of directories and files need to be changed here prior to running any of the notebooks in the analysis/ folder

#### BULK DATASETS ####
# Raw files - Change to location of files after downloading from respective databases
bseq = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq")
bspan = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Kang/genes_matrix_csv")
gtex = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx")
psych = file.path("/home/neuro/Documents/BrainData/Bulk/PsychEncode")

# Outfiles

### Change path to outfile in the file.path() variable
bseq.out = file.path("/home/neuro/Documents/BrainData/Bulk/Brainseq/Formatted")
bspan.out = file.path("/home/neuro/Documents/BrainData/Bulk/BrainSpan/Formatted")
gtex.out = file.path("/home/neuro/Documents/BrainData/Bulk/GTEx/Formatted")
pe.out = file.path("/home/neuro/Documents/BrainData/Bulk/PsychEncode/Formatted")
hdbr.out = file.path("/home/neuro/Documents/BrainData/Bulk/HDBR/Formatted")



#### Single-nucleus datasets ####

# Raw files
aldinger = file.path("/home/neuro/Documents/BrainData/single-cell/aldinger")
hca = file.path("/home/neuro/Documents/BrainData/single-cell/hca/Raw")

# Outfiles

aldinger.out = file.path("/home/neuro/Documents/BrainData/single-cell/aldinger")
