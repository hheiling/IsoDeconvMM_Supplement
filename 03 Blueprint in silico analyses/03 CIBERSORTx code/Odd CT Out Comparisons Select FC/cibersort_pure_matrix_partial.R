# Create pure count matrices for CIBERSORTx input
## Goal of analysis: Compare accuracy of CIBERSORTx method when different numbers of genes are used 
## Using top cluster results from 'cibersort_DESeq2.R', subset pure counts matrix from 
## 'cibersort_pure_matrix.R' with 100, 50, and 25 clusters (same subset of clusters used to 
## create mixture matrices in 'cibersort_mix_matrix.R')

library(stringr)

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
# header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Simulation results prefix
prefix_ciber = str_c(header_nas,"CIBERSORTx")

# load clusters with differential expression across cellTypes
# List from Odd CT Out Comparisons Select FC / cibersort_DESeq2.R: 
# Includes the top 100 clusters from DESeq2 output
# as well as subsets of 50 and 25 used to create mixture matrix subsets
# (Top clusters: clusters upregulated in one cell type compared to the other two cell types combined)
# Object: best_clusters list object
load(sprintf("%s/Best_DiffExp_ClusterSubsets_Upregulated.RData", prefix_ciber))
clust_list = best_clusters

# Read in full pure counts matrix from 'cibersort_pure_matrix.R' output
X = read.table(sprintf("%s/Pure_Count_Matrix.txt", prefix_ciber))

# Subset with top 100 clusters from DESeq2 output
r = which(X[,1] %in% clust_list[["100"]])
X_100 = X[c(1,r),]

# Subset of 50 from top 100
r = which(X[,1] %in% clust_list[["50"]])
X_50 = X[c(1,r),]

# Subset of 25 from top 100
r = which(X[,1] %in% clust_list[["25"]])
X_25 = X[c(1,r),]

# Subset of 10 from top 100
r = which(X[,1] %in% clust_list[["10"]])
X_10 = X[c(1,r),]

write.table(X_100, file = sprintf("%s/Pure_Count_Matrix_100.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_50, file = sprintf("%s/Pure_Count_Matrix_50.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_25, file = sprintf("%s/Pure_Count_Matrix_25.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_10, file = sprintf("%s/Pure_Count_Matrix_10.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
# Note: Need quote = F and sep = "\t" in the write.table line for CIBERSORTx to read data properly.

#####################################################################################################

# Create phenotype file
## Column 1 = cell types (correspond to column names of pure count matrices)
## Other columns correspond to each sample present in the pure count matrix (also columns)
## Entries: 1 if column corresponds to cell type, 2 if wanting to compare sample with cell 
## type (do not correspond to given cell type)
pheno_mat = matrix(0, nrow = 3, ncol = 15+1)
pheno_mat[,1] = c("CT1","CT2","CT3")

ct1_vec = c(1,2,2)
ct2_vec = c(2,1,2)
ct3_vec = c(2,2,1)

ct1_block = matrix(ct1_vec, nrow = 3, ncol = 5, byrow = F)
ct2_block = matrix(ct2_vec, nrow = 3, ncol = 5, byrow = F)
ct3_block = matrix(ct3_vec, nrow = 3, ncol = 5, byrow = F)

pheno_mat[,-1] = cbind(ct1_block, ct2_block, ct3_block)

noquote(pheno_mat)

write.table(pheno_mat, file = sprintf("%s/Phenotype_Classes_Blueprint.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")

#####################################################################################################

print(gc())

q(save = 'no')

#####################################################################################################