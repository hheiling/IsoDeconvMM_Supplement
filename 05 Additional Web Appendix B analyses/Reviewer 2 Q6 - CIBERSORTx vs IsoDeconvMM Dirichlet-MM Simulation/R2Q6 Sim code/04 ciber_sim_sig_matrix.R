
# Using the pure cell type count matrices created using 'ciber_sim_pure_matrix.R', 
# create signature matrices to use in CIBERSORTx.

library(stringr)

path_orig = "~/Longleaf/deconvolution/CIBERSORTx_Sim/count matrices/"
path_new = "~/Longleaf/deconvolution/CIBERSORTx_Sim/signature matrices"
files = list.files(path = path_orig, pattern = "Pure_Count_Matrix", full.names = T)

for(f in 1:length(files)){
  # Extract pure cell type count matrix
  X = read.table(files[f])
  x_vec = as.numeric(as.matrix(X[-1,-1]))
  x_mat = matrix(x_vec, nrow = nrow(X)-1, ncol = ncol(X)-1)
  
  # Normalize the counts
  ## average count total for the samples
  total_cts = mean(colMeans(x_mat))
  ## divide by sample count total and multiply by total_cts (no rounding)
  x_norm = x_mat / matrix(colMeans(x_mat), nrow = nrow(x_mat), ncol = ncol(x_mat), byrow = T) * total_cts
  
  # Take the average of the normalized gene counts across the samples (round to nearest integer)
  x_sig = matrix(0, nrow = nrow(x_norm), ncol = 3)
  for(ct in 1:3){
    if(ct == 1){
      idx = 1:5
    } 
    if(ct == 2){
      idx = 6:10
    } 
    if(ct == 3){
      idx = 11:15
    }
    
    x_sig[,ct] = round(rowMeans(x_norm[,idx]))
  }
  
  x_sig = cbind(as.character(X[-1,1]), x_sig)
  x_sig = rbind(c("clusts","CT1","CT2","CT3"),x_sig)
  print(head(x_sig))
  
  file_name = str_remove(basename(files[f]), ".txt")
  file_name = str_c(file_name,"_CiberSig.txt")

  write.table(x_sig, file = sprintf("%s/%s", path_new, file_name),
              row.names = F, col.names = F, quote = F, sep = "\t")
}

############################################################################################
# Experimentation
############################################################################################


X = read.table("~/Longleaf/Blueprint/CIBERSORTx_OddCT_UpregFC/Pure_Count_Matrix_100.txt")


X_new = read.table(sprintf("%s/Pure_Count_Matrix_100_CiberSim.txt",path_new))

X_Chong = read.table(sprintf("%s/Chong_Examples/CIBERSORTx_input_observed_TPM.txt",path_new))

X_ChongSig = read.table(sprintf("%s/Chong_Examples/CIBERSORTx_input_signature_gene_TPM.txt",path_new))

############################################################################################
# X_pre = read.table(sprintf("%s/Pure_Count_Matrix_100.txt",path_orig))
# X_post = X_pre
# num = c(str_c("00",1:9), str_c("0",10:99),"100")
# X_post[,1] = c("geneSymbol",str_c("gene",num))
# 
# write.table(X_post, file = sprintf("%s/Pure_Count_Matrix_100_Test.txt", path_new),
#             row.names = F, col.names = F, quote = F, sep = "\t")
# 
# Sig100_test = X_pre[,c(1,2,7,12)]
# write.table(Sig100_test, file = sprintf("%s/Sig100_Test.txt", path_new),
#             row.names = F, col.names = F, quote = F, sep = "\t")
############################################################################################

############################################################################################


load("~/Longleaf/deconvolution/Simulated_Output/Ciber_Sim_MixProbCombos100.RData")

############################################################################################

############################################################################################

# for(f in 1:length(files)){
#   X = read.table(files[f])
#   # print(any(is.na(X)))
#   # print(any(apply(X[-1,-1], MARGIN = 2, var) < 0.01))
#   # print(any(apply(X[-1,-1], MARGIN = 1, var) < 0.01))
#   # print(summary(apply(X[-1,-1], MARGIN = 2, sd)))
#   # print(summary(apply(X[-1,-1], MARGIN = 1, sd)))
#   # print(length(unique(X[-1,1])))
#   # svd_test = svd(scale(as.matrix(X[-1,-1])))
#   x_vec = as.numeric(as.matrix(X[-1,-1]))
#   x_mat = matrix(x_vec, nrow = nrow(X)-1, ncol = ncol(X)-1)
#   x_scale = scale(x_mat)
#   test_svd = svd(x_scale)
#   
#   X_new = cbind(as.character(X[-1,1]), x_mat)
#   print(head(X_new))
#   X_new = rbind(X[1,], X_new)
#   print(head(X_new))
#   
#   file_name = str_remove(basename(files[f]), ".txt")
#   file_name = str_c(file_name,"_CiberSim.txt")
#   
#   write.table(X_new, file = sprintf("%s/%s", path_new, file_name),
#               row.names = F, col.names = F, quote = F, sep = "\t")
# }