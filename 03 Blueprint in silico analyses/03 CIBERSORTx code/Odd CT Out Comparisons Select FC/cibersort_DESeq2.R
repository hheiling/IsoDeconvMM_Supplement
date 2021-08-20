# Run DESeq2 on pure samples count matrix created in 'cibersort_pure_matrix.R'
# Goal: Find differentially expressed clusters for use in CIBERSORTx

library(DESeq2)
library(stringr)

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

prefix_ciber = str_c(header_nas,"CIBERSORTx")

dat = read.table(sprintf("%s/Pure_Count_Matrix.txt", prefix_ciber), header = T, sep = " ", as.is = T)

clust_names = dat[-1,1]
cellTypes = dat[1,-1]
counts_mat = matrix(as.numeric(as.matrix(dat[-1,-1])), nrow = (nrow(dat)-1), ncol = (ncol(dat)-1))
rownames(counts_mat) = clust_names
colnames(counts_mat) = cellTypes

group = c(rep(1,5), rep(0,10))
# Note: first column = names of clusters
columns1v23 = 1:15
columns2v13 = c(6:10,1:5,11:15)
columns3v12 = c(11:15,1:10)

#Function
comp = function(data,columns,group,title_table){
  
  countmatrix = round(data[,columns])
  colheaders = as.matrix(colnames(countmatrix))
  colinfo = cbind(colheaders,group)
  dds = DESeqDataSetFromMatrix(countData = countmatrix, colData = colinfo, design = ~ group)
  
  dds_analysis = DESeq(dds)
  # results = results(dds_analysis, alpha = 0.05)
  results = results(dds_analysis)
  
  siggenes = results[order(results$padj),]
  write.table(siggenes, file = title_table, col.names = NA, row.names = T, quote = F, sep = "\t")
  
  return(siggenes)
}

comp1 = comp(data = counts_mat, columns1v23, group, 
               title_table = sprintf("%s/Cibersort_DiffExp_Clusters_CT1_vs_CT23.txt", prefix_ciber))
comp2 = comp(data = counts_mat, columns2v13, group, 
               title_table = sprintf("%s/Cibersort_DiffExp_Clusters_CT2_vs_CT13.txt", prefix_ciber))
comp3 = comp(data = counts_mat, columns3v12, group, 
               title_table = sprintf("%s/Cibersort_DiffExp_Clusters_CT3_vs_CT12.txt", prefix_ciber))


head(comp1)

# Only examine clusters where the log2 fold change is positive
## Note: want clusters that up upregulated in the odd cell type out: group = 1
## Given fold change results are treated (group = 1) / untreated (group = 0)

comp1B = comp1[which(comp1[,2] > 0),]
comp2B = comp2[which(comp2[,2] > 0),]
comp3B = comp3[which(comp3[,2] > 0),]

# Find 100 differentially expressed clusters that will differentiate between all three groups

head(comp1B)

n = 34

top_clusts = union(union(rownames(comp1B[1:n,]), rownames(comp2B[1:n,])), 
                       rownames(comp3B[1:n,]))

length(top_clusts) # Find an n where this results in 100 top clusters

set.seed(3658)
top_clusts = sample(top_clusts, 100, replace = F)
save(top_clusts, file = sprintf("%s/Top_DiffExp_Clusters_Upregulated.RData", prefix_ciber))

# Find subsets of 50 and 25 of 'best' upregulated clusters

top_100 = top_clusts
set.seed(2020)
seeds = sample(1000:9999, 10, replace = F)

## Best 50 genes
n = 17
top_50 = union(union(rownames(comp1B[1:n,]), rownames(comp2B[1:n,])), 
               rownames(comp3B[1:n,]))
length(top_50)
set.seed(seeds[1])
top_50 = sample(top_50, 50, replace = F)
length(top_50)

## Best 25 genes
n = 9
top_25 = union(union(rownames(comp1B[1:n,]), rownames(comp2B[1:n,])), 
               rownames(comp3B[1:n,]))
length(top_25)
set.seed(seeds[2])
top_25 = sample(top_25, 25, replace = F)
length(top_25)

## Best 10 genes
n = 4
top_10 = union(union(rownames(comp1B[1:n,]), rownames(comp2B[1:n,])), 
               rownames(comp3B[1:n,]))
length(top_10)
set.seed(seeds[3])
top_10 = sample(top_10, 10, replace = F)
length(top_10)

best_clusters = list()
best_clusters[["100"]] = top_100
best_clusters[["50"]] = top_50
best_clusters[["25"]] = top_25
best_clusters[["10"]] = top_10 

save(best_clusters, file = sprintf("%s/Best_DiffExp_ClusterSubsets_Upregulated.RData",prefix_ciber))

################################################################################################

# Example of running DESeq2 from Wei Sun:
# 
# colData = meta_ind
# 
# for(i in 1:ncol(colData)){
#   
#   if(is.character(colData[[i]])){
#     colData[[i]] = as.factor(colData[[i]])
#   }
#   
# }
# 
# dim(colData)
# 
# colData[1:2,]
# 
# summary(colData)
# 
# colData$diagnosis = factor(colData$diagnosis, levels=c("Control", "ASD"))
# 
# dd0 = DESeqDataSetFromMatrix(countData = trec1,
#                              colData = colData,
#                              design = ~ diagnosis)
# 
# dd0  = DESeq(dd0)
# 
# res0 = results(dd0)
# 
# dim(res0)
# 
# res0[1:2,]