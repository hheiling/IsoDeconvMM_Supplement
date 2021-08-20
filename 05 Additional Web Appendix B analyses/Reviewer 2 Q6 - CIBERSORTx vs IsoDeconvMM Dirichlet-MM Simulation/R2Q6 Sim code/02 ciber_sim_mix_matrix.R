# Create mixture count matrix for CIBERSORTx input

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

# load genes_df data.frame object
load(sprintf("%s/genes_cibersort.RData", prefix_mat))

# Identify clusters with differential expression
clusts = as.character(genes_df$diffExp)
print(clusts)
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
    r = which(str_detect(exon_cts[,2], clusts[g]))
    # Find total count for gene
    ct = sum(exon_cts[r,1])
    X[g,i] = ct
  }
  
  cat("End file ", basename(files[i]), "\n")
}

# First column: gene names
# First row: names of mixture files
X = cbind(rownames(X), X)
X = rbind(c("clusts",mix_names), X)
## Note: move mix file 100 to very end so it is in the correct order
head(X)
X[1,11]

X = cbind(X,X[,11]) 
X = X[,-11]
rownames(X) = NULL
colnames(X) = NULL

write.table(X, file = sprintf("%s/Mix_Count_Matrix_100.txt", prefix_ciber),
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

# Same seed as in 'ciber_sim_pure_matrix.R'
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

write.table(X_50, file = sprintf("%s/Mix_Count_Matrix_50.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_25, file = sprintf("%s/Mix_Count_Matrix_25.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(X_10, file = sprintf("%s/Mix_Count_Matrix_10.txt", prefix_ciber),
            row.names = F, col.names = F, quote = F, sep = "\t")



#####################################################################################################

print(gc())

q(save = 'no')

#####################################################################################################