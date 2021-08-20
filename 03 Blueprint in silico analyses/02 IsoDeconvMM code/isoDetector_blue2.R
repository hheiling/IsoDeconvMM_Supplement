# Running isoDetector for each pure reference sample

# Changes from original isoDetector_blue.R:
# Changed file split between pure and mixture files
#   Mixture file creation: all pure cell type files must come from same person
# Will use the combined fragment length distribution file instead of individual 
#   (possibly randomly selected) fragment length distribution files

library(stringr)
library(isoform)

# Arrays 1-30
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
batch <- array_val

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Materials prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Pure cell type files prefix (to be used for algorithm fit)
prefix_pure = str_c(header_pine,"Fit_Samples2")

# Datasets for cell types
cellTypes = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")
## Assume ordering of datasets in CT correspond to cell types CT1, CT2, and CT3

# Find pure reference files to use in algorithm fit
CT1_files = list.files(path = str_c(prefix_pure, cellTypes[1], sep = "/"), full.names = T)
CT2_files = list.files(path = str_c(prefix_pure, cellTypes[2], sep = "/"), full.names = T)
CT3_files = list.files(path = str_c(prefix_pure, cellTypes[3], sep = "/"), full.names = T)

set.seed(8500)
CT_ref_files = c(sample(CT1_files, size = 20, replace = F),
                 sample(CT2_files, size = 20, replace = F),
                 sample(CT3_files, size = 20, replace = F))

# Find pure reference samples to use for finding clusters (separate from samples 
# to use in algorithm fit)
# Chose cell type (dataset)
set.seed(8766)
seeds = sample(1000:9999, 3, replace = F)
if(batch <= 10){
  CT = "EGAD00001002671"
  # Simulated output prefix
  prefix_out = str_c(header_nas,"isoDetector_out2/",CT)
  if(!dir.exists(prefix_out)){dir.create(prefix_out, recursive = T)}
  j = batch
  set.seed(seeds[1])
  pure_files = sample(CT1_files[-which(CT1_files %in% CT_ref_files[1:20])], size = 10, replace = F)
}else if(batch <= 20){ # batch 11 to 20
  CT = "EGAD00001002674"
  # Simulated output prefix
  prefix_out = str_c(header_nas,"isoDetector_out2/",CT)
  if(!dir.exists(prefix_out)){dir.create(prefix_out, recursive = T)}
  j = batch - 10 # want j from 1 to 10
  set.seed(seeds[2])
  pure_files = sample(CT2_files[-which(CT2_files %in% CT_ref_files[21:40])], size = 10, replace = F)
}else if(batch >= 21){ # batch 21 to 30
  CT = "EGAD00001002675"
  # Simulated output prefix
  prefix_out = str_c(header_nas,"isoDetector_out2/",CT)
  if(!dir.exists(prefix_out)){dir.create(prefix_out, recursive = T)}
  j = batch - 20 # want j from 1 to 10 
  set.seed(seeds[3])
  pure_files = sample(CT3_files[-which(CT3_files %in% CT_ref_files[41:60])], size = 10, replace = F)
}

print(basename(pure_files))

pure_select = pure_files[j] # Select specific file to run isoDetector on
print(basename(pure_select))

samp_name = unlist(str_split(basename(pure_select), "_"))[1]
print(samp_name)

# Use fragment length distribution file that is composed of a merging of 50 fragment length distribution files
fragSizeFile = str_c(prefix_mat, "Combined_FragDist.txt", sep="/")

# Find BED and knownIsoforms objects
bedFile = sprintf("%s/gencode.v15.nonoverlap.exon.bed", prefix_mat)
knownIsoforms = sprintf("%s/gencode.v15.nonoverlap.exon.knownIsoforms.RData", prefix_mat)

# Specify name of output file
output_file = str_c(samp_name, "_geneModel_knownIsoforms.RData")

print(output_file)

isoDetector(pure_select, bedFile, fragSizeFile, readLen = 100, # readLen from step1_get_counts.Rout
            sprintf("%s/%s", prefix_out, output_file), 
            knownIsoforms=knownIsoforms)

