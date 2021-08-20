# Create pure count matrix for CIBERSORTx input

library(stringr)

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Blueprint materials output prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Pure cell type files prefix
prefix_pure = str_c(header_pine,"Fit_Samples2")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Mixture_Samples2")

# Simulation results prefix
prefix_results = str_c(header_nas,"CIBERSORTx")
if(!dir.exists(prefix_results)){dir.create(prefix_results, recursive = T)}

# load nTE object
load(sprintf("%s/gencode.v15.nTE.RData", prefix_mat))

# Identify clusts from chromosome 1 through 4
## Note: using clusters instead of genes to keep consistent with IsoDeconvMM set-up
clusts = c()
chrom = c("chr1_","chr2_","chr3_","chr4_")
for(i in 1:length(chrom)){
  clusts = union(clusts, nTE$clustID[which(str_detect(nTE$clustID, chrom[i]))])
}
length(clusts)

# Check that clusts only used once
if(length(clusts) != length(unique(clusts))){
  stop("Make sure clusts are not repeated")
}

# Identify pure samples to use

# Datasets for cell types
CT = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")
## Assume ordering of datasets in CT correspond to cell types CT1, CT2, and CT3
cellTypes = c("CT1","CT2","CT3")

# Find pure reference files
CT1_files = list.files(path = str_c(prefix_pure, CT[1], sep = "/"), full.names = T)
CT2_files = list.files(path = str_c(prefix_pure, CT[2], sep = "/"), full.names = T)
CT3_files = list.files(path = str_c(prefix_pure, CT[3], sep = "/"), full.names = T)

set.seed(8500) # Seed from Sim 10 series
CT_ref_files = c(sample(CT1_files, size = 20, replace = F),
                 sample(CT2_files, size = 20, replace = F),
                 sample(CT3_files, size = 20, replace = F))

pure_ref_files_all = cbind(CT_ref_files, c(rep("CT1",times=20),rep("CT2",times=20),rep("CT3",times=20)))

# Select first 5 pure samples from each cell type
## Note: Same samples used in Sim 10A
pure_ref_files = pure_ref_files_all[c(1:5,21:25,41:45),]

cat("dim pure_ref_files: ",dim(pure_ref_files),"\n")

files = pure_ref_files[,1]
cellTypes = pure_ref_files[,2]

# Initialize count matrix
X = matrix(0, nrow = length(clusts), ncol = length(files))
rownames(X) = clusts
colnames(X) = cellTypes # CT1, CT2, or CT3

for(i in 1:length(files)){
  # Read in sample
  ## Single column of counts, rownames = exon sets
  exon_cts = read.table(files[i], as.is = T)
  
  for(g in 1:length(clusts)){
    # Identify rows that correspond to gene of interest
    r = which(str_detect(rownames(exon_cts), clusts[g]))
    # Find total count for gene
    ct = sum(exon_cts[r,])
    X[g,i] = ct
  }
  
  cat("End file ", basename(files[i]), "\n")
}

# First column: cluster names
# First row: names of cell types
X = cbind(rownames(X), X)
X = rbind(c("clusts", pure_ref_files[,2]), X)
rownames(X) = NULL
colnames(X) = NULL

write.table(X, file = sprintf("%s/Pure_Count_Matrix.txt", prefix_results))

#####################################################################################################

print(gc())

q(save = 'no')

#####################################################################################################