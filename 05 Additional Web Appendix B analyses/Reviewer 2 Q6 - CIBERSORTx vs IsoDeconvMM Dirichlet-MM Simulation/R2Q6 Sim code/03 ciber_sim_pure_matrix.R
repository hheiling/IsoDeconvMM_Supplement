# Create pure count matrix for CIBERSORTx input

library(stringr)

header_nas = "/nas/longleaf/home/hheiling/"
header_pine = "/pine/scr/h/h/hheiling/"

# Prefix for needed materials
prefix_mat = str_c(header_pine,"Simulated_Output")
# Pure cell type files prefix
prefix_pure = str_c(header_pine,"Simulated_Files/Ciber_Sim_Pure_Sample_Counts")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header_pine,"Simulated_Files/Ciber_Sim_Mixture_Sample_Counts")

# Simulation results prefix
prefix_ciber = str_c(header_nas,"CIBERSORTx_Sim")
if(!dir.exists(prefix_ciber)){dir.create(prefix_ciber, recursive = T)}

# load genes_df data.frame object
load(sprintf("%s/genes_cibersort.RData", prefix_mat))

# Identify clusters with differential expression
clusts = as.character(genes_df$diffExp)
print(clusts)
length(clusts)

# Identify pure samples to use

# Find pure reference files
CT1_files = list.files(path = prefix_pure, pattern = "CT1", full.names = T)
CT2_files = list.files(path = prefix_pure, pattern = "CT2", full.names = T)
CT3_files = list.files(path = prefix_pure, pattern = "CT3", full.names = T)

# Use the files labels 1-5 of the 20 total files
## Note: files 11-20 used in mixture sample creation
## Only using 5 of the pure reference files per cell type to match Blueprint simultaions
CT_ref_files = c(CT1_files[1:5], CT2_files[1:5], CT3_files[1:5])

pure_ref_files = cbind(CT_ref_files, c(rep("CT1",times=5),rep("CT2",times=5),rep("CT3",times=5)))

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
    r = which(str_detect(exon_cts[,2], clusts[g]))
    # Find total count for gene
    ct = sum(exon_cts[r,1])
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

write.table(X, file = sprintf("%s/Pure_Count_Matrix_100.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")

#####################################################################################################
# Select 50, 25, and 10 genes
#####################################################################################################

# Randomly select approx equal number of genes from the genes up-regulated in CT1, CT2, and CT3
## 50 total genes: select (16,17,17)
## 25 total genes: select (8,8,9)
## 10 total genes: select (3,3,4)
## Make sure the 25 genes are a subest of the 50 genes, and the 10 genes are a subset of the 25 genes

## 33 genes up-regulated in CT1 compared to other cell types
CT1_upreg = as.character(genes_df$diffExp[which(genes_df$CT_UpReg == "CT1")])
print(length(CT1_upreg))
## 33 genes up-regulated in CT2 compared to other cell types
CT2_upreg = as.character(genes_df$diffExp[which(genes_df$CT_UpReg == "CT2")])
print(length(CT2_upreg))
## 34 genes up-regulated in CT3 compared to other cell types
CT3_upreg = as.character(genes_df$diffExp[which(genes_df$CT_UpReg == "CT3")])
print(length(CT3_upreg))

set.seed(2021)
genes_50 = c(sample(CT1_upreg, 16, replace = F),
             sample(CT2_upreg, 17, replace = F),
             sample(CT3_upreg, 17, replace = F))
genes_50 = cbind(genes_50, c(rep("CT1",16),rep("CT2",17),rep("CT3",17)))
print(genes_50)

genes_25 = c(sample(genes_50[which(genes_50[,2] == "CT1"),1], 8),
             sample(genes_50[which(genes_50[,2] == "CT2"),1], 8),
             sample(genes_50[which(genes_50[,2] == "CT3"),1], 9))
genes_25 = cbind(genes_25,c(rep("CT1",8),rep("CT2",8),rep("CT3",9)))
print(genes_25)

genes_10 = c(sample(genes_25[which(genes_25[,2] == "CT1"),1], 3),
             sample(genes_25[which(genes_25[,2] == "CT2"),1], 3),
             sample(genes_25[which(genes_25[,2] == "CT3"),1], 4))
genes_10 = cbind(genes_10, c(rep("CT1",3),rep("CT2",3),rep("CT3",4)))
print(genes_10)


X_50 = X[c(1,which(X[,1] %in% genes_50[,1])),]
X_25 = X[c(1,which(X[,1] %in% genes_25[,1])),]
X_10 = X[c(1,which(X[,1] %in% genes_10[,1])),]

write.table(X_50, file = sprintf("%s/Pure_Count_Matrix_50.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_25, file = sprintf("%s/Pure_Count_Matrix_25.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_10, file = sprintf("%s/Pure_Count_Matrix_10.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")

# Note:
# We can use the same "Phenotype_Classes_Blueprint.txt" created in 'cibersort_pure_matrix_partial.R' 
# code file used for the Blueprint CIBERSORTx analyses

#####################################################################################################

print(gc())

q(save = 'no')

#####################################################################################################