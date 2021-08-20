# Outlier removal processs

X_mat = function(pList){ # clust_names
  
  # q = number genes
  # J = number cell types
  # n = number samples
  
  Xt = pList[[1]]
  # Xt = Xt[which(rownames(Xt) %in% clust_names),]
  
  for(i in 2:length(pList)){
    Xt = cbind(Xt, pList[[i]])
  }
  
  # Dim of Xt: q x (J*n)
  # Dim of t(Xt): (J*n) x q
  X = t(Xt)
  
  return(X)
  
}

outlier_removal = function(X){
  
  # q = number genes
  
  # Distances based on correlation
  ## qxq matrix where all diagonals = 0
  D = as.matrix(as.dist(1-(cor(X)+1)/2))
  
  # Find median distances for each gene
  M = numeric(nrow(D))
  
  for(i in 1:nrow(D)){
    # Leave out diagonal = 0
    M[i] = median(D[i,-i])
  }
  
  names(M) = rownames(D)
  
  q25 = quantile(M, probs = 0.25)
  q75 = quantile(M, probs = 0.75)
  IQR = q75 - q25
  cutoff = q75 + 1.5*IQR
  
  frac_out = mean(M > cutoff)
  
  M_order = M[order(M, decreasing = T)]
  if(frac_out < 0.10){
    clust_keep = names(M_order[-(1:round(0.1*length(M_order)))])
  }else if(frac_out > 0.20){
    clust_keep = names(M_order[-(1:round(0.2*length(M_order)))])
  }else{
    clust_keep = names(M[which(M > cutoff)])
  }
  
  return(clust_keep)
  
}

outlier_clustering = function(X, cutree_method = c("k","h")){
  
  # q = number genes
  q = ncol(X)
  
  # Distances based on correlation
  D = as.dist(1-(cor(X)+1)/2)
  
  H = hclust(D)
  
  if(length(cutree_method) > 1){
    method = cutree_method[1]
  }else{
    method = cutree_method
  }
  
  if(!(method %in% c("k","h"))){
    stop("cutree_method must be either 'k' or 'h'")
  }
  
  if(method == "h"){
    # Option A
    h_seq = seq(0.2, 0.8, by=0.05)

    C_lst = list()
    prop_main = numeric(length(h_seq))

    for(h in 1:length(h_seq)){
      C = cutree(H, h = h_seq[h])
      C_lst[[h]] = C

      # Find size of the largest cluster for specified h
      prop_main[h] = max(table(C) / q)
    }

    # Ignore cases when only a single cluster produced
    if(any(prop_main == 1)){
      prop_main[which(prop_main == 1)] = 0
    }

    if(max(prop_main) < 0.75){
      C_keep = C_lst[[which.max(prop_main)]]
    }else{
      options = which(prop_main >= 0.75)
      C_keep = C_lst[[min(options)]]
    }

    # Find the clust_names in the largest cluster of C_keep
    k_nums = table(C_keep)
    k_keep = as.numeric(names(k_nums[which.max(k_nums)]))
    clust_names = names(C_keep[which(C_keep == k_keep)])
  }else if(method == "k"){
    # Option B
    C = cutree(H, k=2)
    prop_k = table(C) / q
    k_keep = which.max(prop_k)
    clust_names = names(C[which(C == k_keep)])
  }
  
  
  return(clust_names)
  
}

################################################################################################

p_dfA = function(pList, mix_labels, clust_names = NULL, p_combos, out_remove = F){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  if(is.null(clust_names)){
    clust_names = rownames(pList[[1]])
  }
  
  if(out_remove == T){
    # Find (J*n) x q matrix of prob estimates for all samples and cell types for the q genes
    X_full = X_mat(pList = pList)
    
    # Restrict X to only include genes from clusters of interest
    X = X_full[,which(colnames(X_full) %in% clust_names)]
    
    
    # Based on correlation distances, remove outliers and output new cluster restrictions
    cat("length of original clust_names: ", length(clust_names), "\n")
    clust_names = outlier_removal(X = X)
    cat("length of new clust_names after filtering method 1: ", length(clust_names), "\n")
  }
  
  # Calculate overall prob estimates for each sample based on gene-level prob estimates 
  # from clusters of interest
  p_est = matrix(NA, nrow = length(pList), ncol = J)
  
  for(i in 1:length(pList)){
    p_mat = pList[[i]]
    p.fin = spatial.median(X = p_mat[which(p_mat$WARN <= 1 & rownames(p_mat) %in% clust_names),-c(J,(J+1))])
    fin_est = c(p.fin,1-sum(p.fin))
    
    p_est[i,] = fin_est
  }
  
  colnames(p_est) = c("CT1","CT2","CT3")
  rownames(p_est) = mix_labels
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}

################################################################################################

# Remove clusters with extreme probabilities (near 0 and 1 for all cell type proportion estimates)
p_dfA_out2 = function(pList, mix_labels, clust_names = NULL, p_combos){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  if(is.null(clust_names)){
    clust_names = rownames(pList[[1]])
  }
  
  # Calculate overall prob estimates for each sample based on gene-level prob estimates 
  # from clusters of interest
  p_est = matrix(NA, nrow = length(pList), ncol = J)
  
  num_clusts_keep = numeric(length(pList))
  
  for(i in 1:length(pList)){
    p_mat_full = pList[[i]]
    
    # Remove very extreme proportion results
    extreme = ifelse(p_mat_full > 0.99, 1, 0)
    row_keep = which(rowSums(extreme) == 0)
    num_clusts_keep[i] = length(row_keep)
    p_mat = p_mat_full[row_keep,]
    
    p.fin = spatial.median(X = p_mat[which((p_mat$WARN <= 1) & (rownames(p_mat) %in% clust_names)),-c(J,(J+1))])
    fin_est = c(p.fin,1-sum(p.fin))
    
    p_est[i,] = fin_est
  }
  
  cat("Summary of clusters kept after filtering out extreme estimates: \n")
  print(summary(num_clusts_keep))
  
  colnames(p_est) = c("CT1","CT2","CT3")
  rownames(p_est) = mix_labels
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}

################################################################################################

p_dfA_out3 = function(pList, mix_labels, clust_names = NULL, p_combos, out_clust = T,
                      cutree_method = c("k","h")){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  if(is.null(clust_names)){
    clust_names = rownames(pList[[1]])
  }
  
  if(out_clust == T){
    # Find (J*n) x q matrix of prob estimates for all samples and cell types for the q genes
    X_full = X_mat(pList = pList)
    
    # Restrict X to only include genes from clusters of interest
    X = X_full[,which(colnames(X_full) %in% clust_names)]
    
    # Based on correlation distances, remove outliers and output new cluster restrictions
    cat("length of original clust_names: ", length(clust_names), "\n")
    clust_names = outlier_clustering(X = X, cutree_method = cutree_method)
    cat("length of new clust_names after filtering method 3: ", length(clust_names), "\n")
  }
  
  # Calculate overall prob estimates for each sample based on gene-level prob estimates 
  # from clusters of interest
  p_est = matrix(NA, nrow = length(pList), ncol = J)
  
  for(i in 1:length(pList)){
    p_mat = pList[[i]]
    p.fin = spatial.median(X = p_mat[which(p_mat$WARN <= 1 & rownames(p_mat) %in% clust_names),-c(J,(J+1))])
    fin_est = c(p.fin,1-sum(p.fin))
    
    p_est[i,] = fin_est
  }
  
  colnames(p_est) = c("CT1","CT2","CT3")
  rownames(p_est) = mix_labels
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}

################################################################################################

p_dfB = function(pList1, pList2, mix_labels, clust_names, p_combos, out_remove = F){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  if(out_remove == T){
    # Find (J*n) x q matrix of prob estimates for all samples and cell types for the q genes
    X1 = X_mat(pList = pList1)
    X2 = X_mat(pList = pList2)
    X_full = cbind(X1, X2)
    
    # Restrict X to only include genes from clusters of interest
    X = X_full[,which(colnames(X_full) %in% clust_names)]
    
    # Based on correlation distances, remove outliers and output new cluster restrictions
    cat("length of original clust_names: ", length(clust_names), "\n")
    clust_names = outlier_removal(X = X)
    cat("length of new clust_names: ", length(clust_names), "\n")
  }
  
  # Calculate overall prob estimates for each sample based on gene-level prob estimates 
  # from clusters of interest. Assume length(pList1) = length(pList2)
  p_est = matrix(NA, nrow = length(pList1), ncol = J)
  
  for(i in 1:length(pList1)){
    p_mat1 = pList1[[i]]
    p_mat2 = pList2[[i]]
    p_mat = rbind(p_mat1, p_mat2)
    p.fin = spatial.median(X = p_mat[which(p_mat$WARN <= 1 & rownames(p_mat) %in% clust_names),-c(J,(J+1))])
    fin_est = c(p.fin,1-sum(p.fin))
    
    p_est[i,] = fin_est
  }
  
  colnames(p_est) = c("CT1","CT2","CT3")
  rownames(p_est) = mix_labels
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}

p_dfC = function(p_results, CT_cols, p_combos){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  p_est = p_results[,CT_cols]
  rownames(p_est) = p_results$Mixture
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}


################################################################################################

# Creating correlation and SSE tables
# Assumptions:
## data.frames have columns CellType, p_est, and True_p

corr_SSE = function(df_list, col_labels, CT){
  
  corr_mat = matrix(0, nrow = length(CT), ncol = length(df_list))
  SSE_mat = matrix(0, nrow = length(CT), ncol = length(df_list))
  
  for(k in 1:length(df_list)){
    df = df_list[[k]]
    for(ct in 1:length(CT)){
      df_ct = df[which(df$CellType == CT[ct]),]
      corr_mat[ct,k] = cor(df_ct$p_est, df_ct$True_p)
      SSE_mat[ct,k] = sum((df_ct$p_est - df_ct$True_p)^2)
    }
  }
  
  colnames(corr_mat) = col_labels
  colnames(SSE_mat) = col_labels
  
  rownames(corr_mat) = CT
  rownames(SSE_mat) = CT
  
  output = list()
  output[["corr"]] = corr_mat
  output[["SSE"]] = SSE_mat
  
  return(output)
  
}

# Converting correlation and SSE tables into data.frames appropriate for ggplot2 bar plots

# Assumptions:
## rownames = cell types
## colnames associated with comparison difference of interest
bar_df = function(table, col_labels, var_name, value_name){
  
  dfA = data.frame(CellType = rownames(table), table)
  rownames(dfA) = NULL
  colnames(dfA) = c("CellType",col_labels)
  
  dfB = melt(dfA, id = "CellType", variable.name = var_name, value.name = value_name)
  
  if(sum(str_detect(dfB[,var_name], "^X")) == nrow(dfB)){
    dfB[,var_name] = str_remove(dfB[,var_name], "^X")
  }
  
  return(dfB)
}

################################################################################################

# Selecting genes based on SSE across all mixture samples
p_dfD1 = function(pList, mix_labels, clust_names = NULL, p_combos,
                  type = c("A","B"), select_num = 75, cut_off = 0.03){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  if(is.null(clust_names)){
    clust_names = rownames(pList[[1]])
  }
  
  # Calculate SSE for each gene in each mixture sample
  SSE = matrix(0, nrow = length(clust_names), ncol = length(pList))
  
  for(i in 1:length(pList)){
    p_mat = pList[[i]]
    # Only consider genes in clust_names
    p_mat = p_mat[which(rownames(p_mat) %in% clust_names),]
    # Calculate SSE for each gene
    for(r in 1:nrow(p_mat)){
      SSE[r,i] = sum((p_mat[r,] - p_combos[i,])^2)
    }
  }
  rownames(SSE) = clust_names
  colnames(SSE) = mix_labels
  
  # Find genes where SSE across samples is lowest or below a certain threshold
  sum_SSE = rowSums(SSE)
  names(sum_SSE) = clust_names
  if(type == "A"){
    rows_keep = which(sum_SSE <= cut_off)
  }else if(type == "B"){
    SSE_mat = cbind(sum_SSE, order(sum_SSE))
    rows_keep = which(SSE_mat[,2] <= select_num)
  }
  
  if(length(rows_keep) <= 1){
    stop("Selected 1 or fewer genes")
  }
  
  clust_names = names(sum_SSE[rows_keep])
  
  if(type == "A"){
    cat("Number clusters used: ", length(clust_names), "\n")
  }
  
  # Calculate overall prob estimates for each sample based on gene-level prob estimates 
  # from clusters of interest
  p_est = matrix(NA, nrow = length(pList), ncol = J)
  
  for(i in 1:length(pList)){
    p_mat = pList[[i]]
    # Only consider genes in clust_names
    p_mat = p_mat[which(rownames(p_mat) %in% clust_names),]
    p.fin = spatial.median(X = p_mat[which(p_mat$WARN <= 1),-c(J,(J+1))])
    fin_est = c(p.fin,1-sum(p.fin))
    
    p_est[i,] = fin_est
  }
  
  colnames(p_est) = c("CT1","CT2","CT3")
  rownames(p_est) = mix_labels
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}

# Selecting genes based on SSE for individual mixture samples at a time
p_dfD2 = function(pList, mix_labels, clust_names = NULL, p_combos,
                  type = c("A","B"), select_num = 75, cut_off = 0.0003){
  
  # q = number genes
  # J = number cell types
  J = 3
  # n = number samples
  
  if(is.null(clust_names)){
    clust_names = rownames(pList[[1]])
  }
  
  # Calculate overall prob estimates for each sample based on gene-level prob estimates 
  # from clusters of interest
  p_est = matrix(NA, nrow = length(pList), ncol = J)
  
  for(i in 1:length(pList)){
    p_mat = pList[[i]]
    # Only consider genes in clust_names
    p_mat = p_mat[which(rownames(p_mat) %in% clust_names),]
    # Calculate SSE for each gene
    SSE = numeric(nrow(p_mat))
    for(r in 1:nrow(p_mat)){
      SSE[r] = sum((p_mat[r,] - p_combos[i,])^2)
    }
    SSE_mat = cbind(SSE, order(SSE))
    
    if(type == "A"){
      p_mat = p_mat[which(SSE_mat[,1]<= cut_off),]
    }else if(type == "B"){
      p_mat = p_mat[which(SSE_mat[,2] <= select_num),]
    }
    
    if(type == "A"){
      cat("Number clusters used: ", nrow(p_mat), "\n")
    }
    
    if(nrow(p_mat) > 1){
      p.fin = spatial.median(X = p_mat[which(p_mat$WARN <= 1),-c(J,(J+1))])
    }else if(nrow(p_mat) == 0){
      cat("No genes selected in mixture ", i, "\n")
      p.fin = rep(1/J, times = (J-1))
    }else{
      cat("Only one gene selected for mixture ", i, "\n")
      p.fin = p_mat
    }
    
    fin_est = c(p.fin,1-sum(p.fin))
    
    p_est[i,] = fin_est
  }
  
  colnames(p_est) = c("CT1","CT2","CT3")
  rownames(p_est) = mix_labels
  
  # Convert p_est matrix into data.frame to use in ggplots
  dfA = data.frame(mixID = rownames(p_est), p_est)
  rownames(dfA) = NULL
  
  dfB = melt(dfA, id = "mixID", variable.name = "CellType", value.name = "p_est")
  
  df_combos = data.frame(pcID = rownames(p_combos), p_combos)
  rownames(df_combos) = NULL
  
  df_truth = melt(df_combos, id = "pcID", variable.name = "CellType", value.name = "True_p")
  
  df = dfB
  df$True_p = df_truth$True_p
  
  return(df)
  
}

################################################################################################
