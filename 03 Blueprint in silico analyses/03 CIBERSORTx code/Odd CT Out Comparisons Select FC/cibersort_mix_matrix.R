# Create mixture count matrix for CIBERSORTx input

library(stringr)

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Blueprint materials output prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Mixture_Samples2")
# Simulation results prefix
prefix_ciber = str_c(header_nas,"CIBERSORTx")

# load nTE object
load(sprintf("%s/gencode.v15.nTE.RData", prefix_mat))

# load clusters with differential expression across cellTypes
# Object: top_clusts (character vector)
load(sprintf("%s/Top_DiffExp_Clusters_Upregulated.RData", prefix_ciber))

clusts = top_clusts
length(clusts)

# Identify mixture samples (100)
files = list.files(path = prefix_mix, full.names = T)
mix_names = str_remove(basename(files),"_counts.txt")

# Initialize count matrix
X = matrix(0, nrow = length(clusts), ncol = length(files))
rownames(X) = clusts
colnames(X) = mix_names

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
  
  
}

# First column: gene names
# First row: names of mixture files
X = cbind(rownames(X), X)
X = rbind(c("clusts",mix_names), X)
rownames(X) = NULL
colnames(X) = NULL

write.table(X, file = sprintf("%s/Mix_Count_Matrix_100.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")


#####################################################################################################

print(gc())

q(save = 'no')

#####################################################################################################