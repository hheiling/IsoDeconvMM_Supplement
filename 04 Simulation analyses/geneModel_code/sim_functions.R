# Functions to simulate gene counts, isoform probability distributions, and exon set counts


#-----------------------------------------------------------------------------#
# Create gene-level output for all J cell types                                #
#-----------------------------------------------------------------------------#

# Create gene_level output for all three cell types
gene_level = function(total_cts, gene_alpha, seed){
  
  set.seed(seed)
  
  # Define variables
  n = nrow(total_cts)
  J = ncol(total_cts)
  CT_names = rep(1:J, each = n)
  ref_num = rep(1:n, times = J)
  
  gene_names = names(gene_alpha)
  
  # Matrix of counts for each gene
  count = matrix(NA, nrow = length(gene_alpha), ncol = n*J)
  
  # Sample from UCLA Dirichlet distribution
  prob = t(rdirichlet(n = n*J, alpha = gene_alpha))
  rownames(prob) = gene_names
  colnames(prob) = str_c("CT",CT_names,":","ref_",ref_num)
  
  # Using above prob vec, sample total_cts from multinomial distribution
  for(j in 1:J){
    for(i in 1:n){
      count[,i+n*(j-1)] = rmultinom(n = 1, size = total_cts[i,j], prob = prob[,i+n*(j-1)])
    }
  }
  
  rownames(count) = rownames(prob)
  colnames(count) = colnames(prob)
  
  # Vector of dispersion parameters
  # For each sample, assign a constant dispersion parameter (applies to all genes)
  ## Parameterization of interest: mu + mu^2 / theta
  theta = sample(90:120, n*J, replace = TRUE)
  names(theta) = colnames(prob)
  
  gene_output = list(p_mat = prob, ct_mat = count, gene_names = names(gene_alpha),
                     theta = theta)
  
  return(gene_output)
} # End gene_level

gene_level2 = function(total_cts, gene_alpha, theta_range, seed){
  
  set.seed(seed)
  
  # Define variables
  n = nrow(total_cts)
  J = ncol(total_cts)
  CT_names = rep(1:J, each = n)
  ref_num = rep(1:n, times = J)
  
  gene_names = names(gene_alpha)
  
  # Matrix of counts for each gene
  count = matrix(NA, nrow = length(gene_alpha), ncol = n*J)
  
  # Sample from UCLA Dirichlet distribution
  prob = t(rdirichlet(n = n*J, alpha = gene_alpha))
  rownames(prob) = gene_names
  colnames(prob) = str_c("CT",CT_names,":","ref_",ref_num)
  
  # Using above prob vec, sample total_cts from multinomial distribution
  for(j in 1:J){
    for(i in 1:n){
      count[,i+n*(j-1)] = rmultinom(n = 1, size = total_cts[i,j], prob = prob[,i+n*(j-1)])
    }
  }
  
  rownames(count) = rownames(prob)
  colnames(count) = colnames(prob)
  
  # Vector of dispersion parameters
  # For each sample, assign a constant dispersion parameter (applies to all genes)
  ## Parameterization of interest: mu + mu^2 / theta
  theta = sample(theta_range, n*J, replace = TRUE)
  names(theta) = colnames(prob)
  
  gene_output = list(p_mat = prob, ct_mat = count, gene_names = names(gene_alpha),
                     theta = theta)
  
  return(gene_output)
}



#-----------------------------------------------------------------------------#
# Select genes for differential expression and differential isoform usage     #
#-----------------------------------------------------------------------------#

diff_genes = function(CT1_counts, nTE_filtered, num_diff = 200, seed){
  
  set.seed(seed)
  
  # Specify vars
  n = ncol(CT1_counts)
  num_genes = nrow(CT1_counts)
  gene_names = rownames(CT1_counts)
  
  # Check inputs
  if(num_diff %% 2 != 0){
    # Note: will split num_diff into diff expression and diff usage. 
    stop("Number of genes selected for differential expression must divisible by 2")
  }
  
  # Subset: 1000 genes of interest from the original 5,172 genes
  counts_subset = CT1_counts[which(gene_names %in% nTE_filtered$geneId),]
  # Restrict genes to those with isoforms between 3 and 15
  genes_nT_limit = nTE_filtered$geneId[which(nTE_filtered$nT <= 15)]
  
  # Find genes of interest from CT1 output with counts above first p25 of counts 
  # (wrt 1000 genes of interest) in at least 90% of the samples 
  q1 = apply(counts_subset, 2, function(x) quantile(x, probs = 0.25))
  above_q1 = matrix(NA, nrow = nrow(counts_subset), ncol = ncol(counts_subset))
  for(j in 1:ncol(counts_subset)){
    above_q1[,j] = ifelse(counts_subset[,j] > q1[j], 1, 0)
  }
  gene_expCut = rownames(counts_subset)[which(rowSums(above_q1) >= n*0.9)]
  
  # Select genes for differential expression
  ## No overlap between diff expression and diff usage genes
  
  # Find intersection of genes with desired number of isoforms and desired expression range
  gene_choices = intersect(gene_expCut,genes_nT_limit)
  cat("number of gene_choices after expression level and isoform number restriction: ",
      length(gene_choices), "\n")
  all_diff = sample(gene_choices, num_diff, replace = F)
  diffExp = sample(all_diff, num_diff/2, replace = F)
  diffUsg = all_diff[-which(all_diff %in% diffExp)]
  noDiff = sample(gene_choices[-which(gene_choices %in% all_diff)], 100, replace = F)
  
  genes_diff = matrix(NA, nrow = num_genes, ncol = 3)
  genes_diff[,1] = ifelse(gene_names %in% diffExp, 1, 0)
  genes_diff[,2] = ifelse(gene_names %in% diffUsg, 1, 0)
  genes_diff[,3] = ifelse(gene_names %in% noDiff, 1, 0)
  
  colnames(genes_diff) = c("diffExp","diffUsg","noDiff")
  rownames(genes_diff) = gene_names
  
  return(genes_diff)
  
}

#-----------------------------------------------------------------------------#
# Apply fold changes to genes for differential expression for specified genes
# from diff_genes() output
#-----------------------------------------------------------------------------#

diff_exp = function(gene_counts, n, J, CT_diffExp = 2, diff_genes_mat, propUp, seed){
  
  set.seed(seed)
  
  # Specify vars
  num_genes = nrow(gene_counts)
  all_genes = rownames(gene_counts)
  
  # Checks
  if(!(CT_diffExp %in% 1:J)){
    stop("CT_diffExp must be a number in 1 to J")
  }else if(length(CT_diffExp) > 1){
    stop("Function ony set up to handle one cell type with diff expression")
  }
  
  # Specify genes selected for differential expression
  diff_idx = which(diff_genes_mat[,"diffExp"] == 1)
  gene_choices = rownames(diff_genes_mat[diff_idx,])
  num_diffExp = sum(diff_genes_mat[,"diffExp"])
  num_upExp = round(num_diffExp*propUp)
  num_downExp = num_diffExp - num_upExp
  
  # Initialize log2 fold change matrix to be all 0
  fc = matrix(matrix(0, nrow = num_genes, ncol = 1))
  
  # Update fold change matrix
  up_diffExp = sample(gene_choices, num_upExp, replace = F)
  down_diffExp = gene_choices[-which(gene_choices %in% up_diffExp)]
  
  fc[which(all_genes %in% up_diffExp),1] = runif(n = num_upExp, min = log2(1.6), max = log2(2))
  fc[which(all_genes %in% down_diffExp),1] = runif(n = num_downExp, min = -log2(2), 
                                                   max = -log2(1.6))
  
  rownames(fc) = all_genes
  colnames(fc) = str_c("CT",CT_diffExp)
  
  # Apply fold change matrix to gene counts of cell type of interest
  
  ## Multiply CT_diffExp counts by 2^(fc + rnorm(1, mean = 0, sd = 0.05))
  gene_cts = gene_counts[,(1+n*(CT_diffExp-1)):(n*CT_diffExp)]
  fc_rand = matrix(fc, nrow = nrow(fc), ncol = n) + matrix(rnorm(n = nrow(fc)*n, mean = 0, sd = 0.05), nrow = nrow(fc))
  gene_cts_fc = gene_cts * 2^(fc_rand)
  ## Determine proportion of counts for gene_choices when no fold change
  propA = colSums(gene_cts[which(all_genes %in% gene_choices),]) / colSums(gene_cts)
  propB = colSums(gene_cts_fc[which(all_genes %in% gene_choices),]) / colSums(gene_cts_fc)
  prop_R = propA / propB
  print("proportion ratio")
  print(prop_R)
  ## Multiply ratio of proportions to the counts affected by fc
  gene_cts_fc[which(fc != 0),] = gene_cts_fc[which(fc != 0),] * 
    matrix(prop_R, nrow = num_diffExp, ncol = n, byrow = T)
  
  gene_counts_new = gene_counts
  gene_counts_new[,(1+n*(CT_diffExp-1)):(n*CT_diffExp)] = round(gene_cts_fc)
  colnames(gene_counts_new) = colnames(gene_counts)
  rownames(gene_counts_new) = rownames(gene_counts)
  
  return(gene_counts_new)
  
}

# diff_exp_ciber = function(gene_counts, n, J, diff_genes_mat, seed,
#                           genes_df){
#   
#   set.seed(seed)
#   
#   # Specify vars
#   num_genes = nrow(gene_counts)
#   all_genes = rownames(gene_counts)
#   
#   # Specify genes selected for differential expression
#   diff_idx = which(diff_genes_mat[,"diffExp"] == 1)
#   gene_choices = rownames(diff_genes_mat[diff_idx,])
#   # num_diffExp = sum(diff_genes_mat[,"diffExp"])
#   
#   gene_counts_new = gene_counts
#   colnames(gene_counts_new) = colnames(gene_counts)
#   rownames(gene_counts_new) = rownames(gene_counts)
#   
#   # Cell type 1:
#   ## Up-regulate the specified 33 genes in CT1 compared to CT2 and CT3
#   
#   for(ct in 1){
#     
#     # Cell type label
#     ct_lab = str_c("CT",ct)
#     # Identify which genes to upreguate in cell type 'ct' 
#     up_diffExp = genes_df$diffExp[which(genes_df$CT_UpReg == ct_lab)]
#     
#     # Initialize log2 fold change matrix to be all 0
#     fc = matrix(matrix(0, nrow = length(up_diffExp), ncol = 1))
#     rownames(fc) = up_diffExp
#     
#     # Fold change values
#     fc[,1] = runif(n = length(up_diffExp), min = log2(1.6), max = log2(2))
#     
#     # Apply fold change matrix to gene counts of cell type of interest
#     
#     ## Multiply CT1 gene counts by 2^(fc + rnorm(1, mean = 0, sd = 0.05))
#     gene_cts = gene_counts[which(all_genes %in% up_diffExp),(1+n*(ct-1)):(n*ct)]
#     fc_rand = matrix(fc, nrow = nrow(fc), ncol = n) + matrix(rnorm(n = nrow(fc)*n, mean = 0, sd = 0.05), nrow = nrow(fc))
#     gene_cts_fc = gene_cts * 2^(fc_rand)
#     
#     gene_counts_new[which(all_genes %in% up_diffExp),(1+n*(ct-1)):(n*ct)] = round(gene_cts_fc)
#     
#   }
#   
#   # Cell types 2 and 3:
#   ## Cell type 2: keep gene expression the same in CT2 for the 33 genes of interest, but
#   ## down-regulate these 33 genes in CT1 and CT3
#   ## Cell type 3: do the same as cell type 2
#   
#   for(ct in 2:3){
#     
#     # Cell type label
#     ct_lab = str_c("CT",ct)
#     # Identify which genes to upreguate in cell type 'ct' 
#     ## These genes will be down-regulated in the other cell types to make sure they are
#     ## up-regulated in cell type 'ct'
#     up_diffExp = genes_df$diffExp[which(genes_df$CT_UpReg == ct_lab)]
#     
#     # Initialize log2 fold change matrix to be all 0
#     fc = matrix(matrix(0, nrow = length(up_diffExp), ncol = 1))
#     rownames(fc) = up_diffExp
#     colnames(fc) = "FC"
#     
#     # Fold change values
#     fc[,1] = runif(n = length(up_diffExp), min = -log2(2), max = -log2(1.6))
#     
#     # Apply fold change matrix to gene counts of cell type of interest
#     
#     ## Multiply gene counts of cell types != ct by 2^(fc + rnorm(1, mean = 0, sd = 0.05))
#     gene_cts = gene_counts[which(all_genes %in% up_diffExp),-c((1+n*(ct-1)):(n*ct))]
#     fc_rand = matrix(fc, nrow = nrow(fc), ncol = n*2) + matrix(rnorm(n = nrow(fc)*n*2, mean = 0, sd = 0.05), nrow = nrow(fc))
#     gene_cts_fc = gene_cts * 2^(fc_rand)
#     
#     gene_counts_new[which(all_genes %in% up_diffExp),-c((1+n*(ct-1)):(n*ct))] = round(gene_cts_fc)
#     
#     
#   }
#   
#   ## Determine proportion of counts for gene_choices when no fold change
#   propA = colSums(gene_counts[which(all_genes %in% gene_choices),]) / colSums(gene_counts)
#   propB = colSums(gene_counts_new[which(all_genes %in% gene_choices),]) / colSums(gene_counts_new)
#   prop_R = propA / propB
#   print("proportion ratio")
#   print(prop_R)
#   
#   return(gene_counts_new)
#   
# }

diff_exp_ciber = function(gene_counts, n, J, diff_genes_mat, seed,
                          genes_df){
  
  set.seed(seed)
  
  # Specify vars
  num_genes = nrow(gene_counts)
  all_genes = rownames(gene_counts)
  
  # Specify genes selected for differential expression
  diff_idx = which(diff_genes_mat[,"diffExp"] == 1)
  gene_choices = rownames(diff_genes_mat[diff_idx,])
  
  gene_counts_new = gene_counts
  colnames(gene_counts_new) = colnames(gene_counts)
  rownames(gene_counts_new) = rownames(gene_counts)
  
  # Cell type ct:
  ## Multiply the gene counts in the other cell types by the constant 0.775
  ## Multiply the gene counts in ct by a fold change value ranging from 1.24 to 1.55
  ## As a result, the fold change between ct and the other cell types is 1.6 to 2.0
  
  down_const = 0.775
  up_fc = c(1.24, 1.55)
  
  for(ct in 1:3){
    
    # Cell type label
    ct_lab = str_c("CT",ct)
    # Identify which genes to upreguate in cell type 'ct' 
    up_diffExp = genes_df$diffExp[which(genes_df$CT_UpReg == ct_lab)]
    
    # Initialize log2 fold change matrix to be all 0
    fc = matrix(matrix(0, nrow = length(up_diffExp), ncol = 1))
    rownames(fc) = up_diffExp
    
    # Fold change values
    fc[,1] = runif(n = length(up_diffExp), min = log2(up_fc[1]), max = log2(up_fc[2]))
    
    # Apply fold change matrix to gene counts of cell type of interest
    
    ## Multiply CT1 gene counts by 2^(fc + rnorm(1, mean = 0, sd = 0.025))
    gene_cts = gene_counts[which(all_genes %in% up_diffExp),(1+n*(ct-1)):(n*ct)]
    fc_rand = matrix(fc, nrow = nrow(fc), ncol = n) + matrix(rnorm(n = nrow(fc)*n, mean = 0, sd = 0.025), nrow = nrow(fc))
    gene_cts_fc = gene_cts * 2^(fc_rand)
    
    gene_counts_new[which(all_genes %in% up_diffExp),(1+n*(ct-1)):(n*ct)] = round(gene_cts_fc)
    
    ## Multiply CT2 and CT3 gene counts by down_const
    down_cts = gene_counts[which(all_genes %in% up_diffExp),-c((1+n*(ct-1)):(n*ct))]
    gene_counts_new[which(all_genes %in% up_diffExp),-c((1+n*(ct-1)):(n*ct))] = round(down_const * down_cts)
    
  }
  
  ## Determine proportion of counts for gene_choices when no fold change
  propA = colSums(gene_counts[which(all_genes %in% gene_choices),]) / colSums(gene_counts)
  propB = colSums(gene_counts_new[which(all_genes %in% gene_choices),]) / colSums(gene_counts_new)
  prop_R = propA / propB
  print("proportion ratio")
  print(prop_R)
  print("summary of proportion ratio")
  print(summary(prop_R))
  
  return(gene_counts_new)
  
}

#-----------------------------------------------------------------------------#
# Determine Ending Fold Change b/w CT ref and CT j                            #
#-----------------------------------------------------------------------------#
calc_diffExp = function(gene_counts_new, gene_counts_orig, diff_genes_mat){
  
  # Define Variables
  
  # Rows associated with genes selected for differential expression
  rows = which(diff_genes_mat[,"diffExp"] == 1)
  # Original gene counts for CT j before fold change applied
  orig = gene_counts_orig[rows,]
  # New gene counts for CT j after fold change and proportion adjustment
  new = gene_counts_new[rows,]
  
  fc_indiv = new / orig
  fc_avg = rowMeans(fc_indiv)
  
  return(list(fc_indiv = fc_indiv, fc_avg = fc_avg))
  
}

#-----------------------------------------------------------------------------#
# Simulate exon set negative binomial means                                   #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms I >= 3 for all genes
# iso_dist = "uniform": all probabilities will be close to 1/I (I = number isoforms)
#     Note: this type only uses max of alphaRange
# iso_dist = "outlier": one probability will be relatively high and the remaining prob
#     will be approx. evenly distributed among the remaining I-1 isoforms
#     Note: this type uses both min and max of alphaRange
# iso_dist = "paired": two probabilities will be relatively high and the remaining
#     probs will be approx. evenly distributed among the remaining I-2 isoforms
#     Note: this type uses both min and max of alphaRange
iso_exon_info = function(genes_info, nTE_filtered, iso_dist, 
                          alphaRange, EffLen_info, 
                          seed = seed){
  
  set.seed(seed)
  
  # names of 1,000 clusters of interest used in simulation
  clust_names = nTE_filtered$clustID
  # names of genes of interest
  ## note: rows of genes_info matrix restricted to genes in nTE_filtered object
  gene_names = rownames(genes_info)
  # Number samples
  n = ncol(genes_info)
  
  # Check
  if(length(gene_names) != length(clust_names)){
    stop("nrow of genes_info does not match nrow of nTE_filtered")
  }
  
  # Check iso_dist
  iso_dist_options = unique(iso_dist)
  if(!all(iso_dist_options %in% c("uniform","outlier","paired"))){
    stop("iso_dist elements must be one of 'uniform', 'outlier', or 'paired'")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_filtered$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # Effective length matrix - ExI (num exon sets x num isoforms)
    X = EffLen_info[[clust]]$X
    # number isoforms
    I = ncol(X)
    # dirichlet alpha parameters for isoforms
    dir_dist = iso_dist[which(names(iso_dist) == gene)]
    if(dir_dist == "uniform"){
      alpha = rep(alphaRange[2], times = I)
    }else if(dir_dist == "outlier"){
      alpha = c(rep(alphaRange[1], times = (I-1)), alphaRange[2])
    }else if(dir_dist == "paired"){
      alpha = c(rep(alphaRange[1], times = (I-2)), rep(alphaRange[2], times = 2))
    }
    
    # isoform probability matrix - Ixn (col = isoform probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    candiIsoform = EffLen_info[[clust]]$candiIsoform
    rownames(rho) = colnames(candiIsoform)
    colnames(rho) = str_c("ref_",1:n)
    
    # scaling factor for gene
    r_g = numeric(n)
    # coefficient for mu_g = X_g %*% beta_g
    beta = matrix(NA, nrow = I, ncol = n)
    
    for(i in 1:n){
      r_g[i] = gene_ct[i] / sum(X %*% rho[,i])
      beta[,i] = rho[,i] * r_g[i]
    }
    
    # Find exon sets corresponding to rows of X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    # negative binomial means for the exon sets within cluster
    # result: each col = neg bin means for sample i of n samples,
    #   each row corresponds with (possible) exon sets
    mu = X %*% beta
    rownames(mu) = exon_sets
    colnames(mu) = str_c("ref_",1:n)
    
    output[[clust]] = list(iso_alpha = alpha, rho = rho, mu = mu, exon_sets = exon_sets)
  }
  
  return(output)
  
}


#-----------------------------------------------------------------------------#
# Simulate exon set negative binomial means - take 2                          #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms I >= 3 for all genes
# iso_dist = "uniform": all probabilities will be close to 1/I (I = number isoforms)
#     Note: this type only uses max of alphaRange
# iso_dist = "outlier": one probability will be relatively high and the remaining prob
#     will be approx. evenly distributed among the remaining I-1 isoforms
#     Note: this type uses both min and max of alphaRange
# iso_dist = "paired": two probabilities will be relatively high and the remaining
#     probs will be approx. evenly distributed among the remaining I-2 isoforms
#     Note: this type uses both min and max of alphaRange
#     Note: In this situation, the isoforms with the highest alpha are different
#       than in the original iso_exon_info() function
iso_exon_info2 = function(genes_info, nTE_filtered, iso_dist, 
                          alphaRange, EffLen_info, 
                          seed = seed){
  
  set.seed(seed)
  
  # names of 1,000 clusters of interest used in simulation
  clust_names = nTE_filtered$clustID
  # names of genes of interest
  ## note: rows of genes_info matrix restricted to genes in nTE_filtered object
  gene_names = rownames(genes_info)
  # Number samples
  n = ncol(genes_info)
  
  # Check
  if(length(gene_names) != length(clust_names)){
    stop("nrow of genes_info does not match nrow of nTE_filtered")
  }
  
  # Check iso_dist
  iso_dist_options = unique(iso_dist)
  if(!all(iso_dist_options %in% c("uniform","outlier","paired","outlier3"))){
    stop("iso_dist elements must be one of 'uniform', 'outlier', or 'paired'")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_filtered$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # Effective length matrix - ExI (num exon sets x num isoforms)
    X = EffLen_info[[clust]]$X
    # number isoforms
    I = ncol(X)
    # dirichlet alpha parameters for isoforms
    dir_dist = iso_dist[which(names(iso_dist) == gene)]
    if(dir_dist == "uniform"){
      alpha = rep(alphaRange[2], times = I)
    }else if(dir_dist == "outlier"){
      alpha = c(rep(alphaRange[1], times = (I-1)), alphaRange[2])
    }else if(dir_dist == "paired"){
      alpha = c(rep(alphaRange[2], times = 2), rep(alphaRange[1], times = (I-2)))
    }
    
    # isoform probability matrix - Ixn (col = isoform probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    candiIsoform = EffLen_info[[clust]]$candiIsoform
    rownames(rho) = colnames(candiIsoform)
    colnames(rho) = str_c("ref_",1:n)
    
    # scaling factor for gene
    r_g = numeric(n)
    # coefficient for mu_g = X_g %*% beta_g
    beta = matrix(NA, nrow = I, ncol = n)
    
    for(i in 1:n){
      r_g[i] = gene_ct[i] / sum(X %*% rho[,i])
      beta[,i] = rho[,i] * r_g[i]
    }
    
    # Find exon sets corresponding to rows of X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    # negative binomial means for the exon sets within cluster
    # result: each col = neg bin means for sample i of n samples,
    #   each row corresponds with (possible) exon sets
    mu = X %*% beta
    rownames(mu) = exon_sets
    colnames(mu) = str_c("ref_",1:n)
    
    output[[clust]] = list(iso_alpha = alpha, rho = rho, mu = mu, exon_sets = exon_sets)
  }
  
  return(output)
  
}


#-----------------------------------------------------------------------------#
# Simulate exon set negative binomial means - take 3                          #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms I >= 3 for all genes
# iso_dist = "uniform": all probabilities will be close to 1/I (I = number isoforms)
#     Note: this type only uses max of alphaRange
# iso_dist = "outlier1": The first isoform (isoform 1 as determined by first column of 
#     knownIsoforms matrix) of the I isoforms will have the highest probability (by a significant 
#     margin) and the remaining isoforms will have small probabilities that are approx. uniform 
#     across these I-1 isoforms. 
#     Note: this type uses both min and max of alphaRange
# iso_dist = "outlier2": The second isoform (isoform 2 as determined by second column of 
#     knownIsoforms matrix) of the I isoforms will have the highest probability (by a significant 
#     margin) and the remaining isoforms will have small probabilities that are approx. uniform 
#     across these I-1 isoforms. 
#     Note: this type uses both min and max of alphaRange
# iso_dist = "outlier2": The third isoform (isoform 3 as determined by third column of 
#     knownIsoforms matrix) of the I isoforms will have the highest probability (by a significant 
#     margin) and the remaining isoforms will have small probabilities that are approx. uniform 
#     across these I-1 isoforms. 
#     Note: this type uses both min and max of alphaRange
iso_exon_info3 = function(genes_info, nTE_filtered, iso_dist, 
                          alphaRange, EffLen_info, 
                          seed = seed){
  
  set.seed(seed)
  
  # names of 1,000 clusters of interest used in simulation
  clust_names = nTE_filtered$clustID
  # names of genes of interest
  ## note: rows of genes_info matrix restricted to genes in nTE_filtered object
  gene_names = rownames(genes_info)
  # Number samples
  n = ncol(genes_info)
  
  # Check
  if(length(gene_names) != length(clust_names)){
    stop("nrow of genes_info does not match nrow of nTE_filtered")
  }
  
  # Check iso_dist
  iso_dist_options = unique(iso_dist)
  if(!all(iso_dist_options %in% c("uniform","outlier1","outlier2","outlier3"))){
    stop("iso_dist elements must be one of 'uniform', 'outlier1', 'outlier2', or 'outlier3'")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_filtered$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # Effective length matrix - ExI (num exon sets x num isoforms)
    X = EffLen_info[[clust]]$X
    # number isoforms
    I = ncol(X)
    # dirichlet alpha parameters for isoforms
    dir_dist = iso_dist[which(names(iso_dist) == gene)]
    if(dir_dist == "uniform"){
      alpha = rep(alphaRange[2], times = I)
    }else if(dir_dist == "outlier1"){
      alpha = c(alphaRange[2], rep(alphaRange[1], times = (I-1)))
    }else if(dir_dist == "outlier2"){
      alpha = c(alphaRange[1], alphaRange[2], rep(alphaRange[1], times = (I-2)))
    }else if(dir_dist == "outlier3"){
      if(I == 3){
        alpha = c(rep(alphaRange[1], times=2), alphaRange[2])
      }else{ # I > 3
        alpha = c(rep(alphaRange[1], times=2), alphaRange[2], rep(alphaRange[1], times = (I-3)))
      }
    }
    
    # isoform probability matrix - Ixn (col = isoform probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    candiIsoform = EffLen_info[[clust]]$candiIsoform
    rownames(rho) = colnames(candiIsoform)
    colnames(rho) = str_c("ref_",1:n)
    
    # scaling factor for gene
    r_g = numeric(n)
    # coefficient for mu_g = X_g %*% beta_g
    beta = matrix(NA, nrow = I, ncol = n)
    
    for(i in 1:n){
      r_g[i] = gene_ct[i] / sum(X %*% rho[,i])
      beta[,i] = rho[,i] * r_g[i]
    }
    
    # Find exon sets corresponding to rows of X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    # negative binomial means for the exon sets within cluster
    # result: each col = neg bin means for sample i of n samples,
    #   each row corresponds with (possible) exon sets
    mu = X %*% beta
    rownames(mu) = exon_sets
    colnames(mu) = str_c("ref_",1:n)
    
    output[[clust]] = list(iso_alpha = alpha, rho = rho, mu = mu, exon_sets = exon_sets)
  }
  
  return(output)
  
}

#-----------------------------------------------------------------------------#
# Simulate exon set negative binomial means - take 4                          #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms 3 <= I <= 15 for all genes
# Alpha parameters taken from Blueprint data results. These parameters must
# be multiplied by 5 (see "Blueprint Alpha Parameter.Rmd" for justification)
iso_exon_info4 = function(genes_info, nTE_filtered, EffLen_info, CT, alpha_record, seed = seed){
  
  set.seed(seed)
  
  # names of 1,000 clusters of interest used in simulation
  clust_names = nTE_filtered$clustID
  # names of genes of interest
  ## note: rows of genes_info matrix restricted to genes in nTE_filtered object
  gene_names = rownames(genes_info)
  # Number samples
  n = ncol(genes_info)
  
  # Check
  if(length(gene_names) != length(clust_names)){
    stop("nrow of genes_info does not match nrow of nTE_filtered")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_filtered$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # Effective length matrix - ExI (num exon sets x num isoforms)
    X = EffLen_info[[clust]]$X
    # number isoforms
    I = ncol(X)
    # Choose alpha from pre-determined alpha_record (already scaled by factor of 5)
    alpha = alpha_record[[clust]][,CT]
    
    # isoform probability matrix - Ixn (col = isoform probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    candiIsoform = EffLen_info[[clust]]$candiIsoform
    rownames(rho) = colnames(candiIsoform)
    colnames(rho) = str_c("ref_",1:n)
    
    # scaling factor for gene
    r_g = numeric(n)
    # coefficient for mu_g = X_g %*% beta_g
    beta = matrix(NA, nrow = I, ncol = n)
    
    for(i in 1:n){
      r_g[i] = gene_ct[i] / sum(X %*% rho[,i])
      beta[,i] = rho[,i] * r_g[i]
    }
    
    # Find exon sets corresponding to rows of X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    # negative binomial means for the exon sets within cluster
    # result: each col = neg bin means for sample i of n samples,
    #   each row corresponds with (possible) exon sets
    mu = X %*% beta
    rownames(mu) = exon_sets
    colnames(mu) = str_c("ref_",1:n)
    
    output[[clust]] = list(iso_alpha = alpha, rho = rho, mu = mu, exon_sets = exon_sets)
  }
  
  return(output)
  
}

#-----------------------------------------------------------------------------#
# Simulate exon set counts assuming Dirichlet-Multinomial Model               #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms 3 <= I <= 15 for all genes
# Alpha parameters taken from Blueprint data results. 
iso_exon_info_MM = function(genes_info, nTE_filtered, EffLen_info, CT, alpha_record, seed = seed){
  
  set.seed(seed)
  
  # names of 1,000 clusters of interest used in simulation
  clust_names = nTE_filtered$clustID
  # names of genes of interest
  ## note: rows of genes_info matrix restricted to genes in nTE_filtered object
  gene_names = rownames(genes_info)
  # Number samples
  n = ncol(genes_info)
  
  # Check
  if(length(gene_names) != length(clust_names)){
    stop("nrow of genes_info does not match nrow of nTE_filtered")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_filtered$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # Effective length matrix - ExI (num exon sets x num isoforms)
    X = EffLen_info[[clust]]$X
    # number isoforms
    I = ncol(X)
    # Choose alpha from pre-determined alpha_record
    ## alpha_record result stores alpha for all three cell types, choose correct column
    alpha = alpha_record[[clust]][,CT]
    
    # isoform probability matrix - Ixn (col = isoform probability vector associated with sample i)
    # let rho = l_tilde * gamma (see Simulation Step Overview to understand notation)
    rho = t(rdirichlet(n = n, alpha = alpha))
    candiIsoform = EffLen_info[[clust]]$candiIsoform
    rownames(rho) = colnames(candiIsoform)
    colnames(rho) = str_c("ref_",1:n)
    
    # Calculate l_tilde so that gamma can be calculated from above rho
    l_tilde = colSums(X)
    gamma_vals = matrix(0, nrow = I, ncol = n) # Same dimensions as rho
    for(i in 1:ncol(rho)){
      gamma_vals[,i] = rho[,i] / l_tilde
    }
    rownames(gamma_vals) = colnames(candiIsoform)
    colnames(gamma_vals) = str_c("ref_",1:n)
    
    # find scaling factor to make r * sum(X %*% gamma) = 1
    r_g = numeric(n)
    
    for(i in 1:n){
      r_g[i] = 1 / sum(X %*% matrix(gamma_vals[,i], ncol = 1))
    }
    
    # Multinomial probabilities: r * X %*% gamma
    mm_probs = matrix(0, nrow = nrow(X), ncol = n)
    for(i in 1:n){
      mm_probs[,i] = r_g[i] * X %*% matrix(gamma_vals[,i], ncol = 1)
    }
    
    # Check
    if(isFALSE(all.equal(colSums(mm_probs), rep(1,times=n)))){
      print(colSums(mm_probs))
      stop("mm_probs do not add up to 1")
    }
    
    # Find exon sets corresponding to rows of X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    rownames(mm_probs) = exon_sets
    colnames(mm_probs) = str_c("ref_",1:n)
    
    output[[clust]] = list(iso_alpha = alpha, rho = rho, gamma_vals = gamma_vals,
                           mm_probs = mm_probs, gene_ct = gene_ct, exon_sets = exon_sets)
  }
  
  return(output)
  
}

#-----------------------------------------------------------------------------#
# Simulate exon set counts for 'other' genes                                  #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms I >= 3 for all genes
# E = number of singular exon sets
# iso_dist = "uniform": all probabilities will be close to 1/E
#     Note: this type only uses max of alphaRange
other_exonset_count = function(genes_info, nTE_other, exon_sets_other, 
                               iso_dist = rep("uniform", times = nrow(nTE_other)), 
                               alphaRange = c(20,50), seed = seed){
  set.seed(seed)
  
  # names of 'other' clusters 
  clust_names = nTE_other$clustID
  # names of 'other' genes
  ## Note: rows of genes_info matrix restricted to genes in nTE_other object
  gene_names = rownames(genes_info)
  # Number samples
  n = ncol(genes_info)
  
  # Check
  if(length(gene_names) != length(clust_names)){
    stop("nrow of genes_info does not match nrow of nTE_other")
  }
  
  # Check iso_dist
  iso_dist_options = unique(iso_dist)
  if(!all(iso_dist_options %in% c("uniform","outlier"))){
    stop("iso_dist elements must be one of 'uniform' or 'outlier'")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_other$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # singular exon sets
    exon_sets = exon_sets_other[[clust]]
    
    # Distribute gene counts to singular exon sets according to iso_dist specification
    ## E = number singular exon sets
    E = length(exon_sets)
    ## dirichlet alpha parameters
    dir_dist = iso_dist[which(gene_names == gene)]
    if(dir_dist == "uniform"){
      alpha = rep(alphaRange[2], times = E)
    }else{
      stop("isoform distribution only allowed to be 'uniform'")
    }
    
    # exon set probability matrix - Exn (col = exon set probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    rownames(rho) = exon_sets
    colnames(rho) = str_c("ref_",1:n)
    
    # Determine exon set counts from multinomial distriution
    exon_set_cts = matrix(NA, nrow = E, ncol = n)
    for(i in 1:n){
      exon_set_cts[,i] = rmultinom(n = 1, size = gene_ct[i], prob = rho[,i])
    }
    
    rownames(exon_set_cts) = exon_sets
    colnames(exon_set_cts) = colnames(genes_info)
    
    output[[clust]] = list(exon_sets = exon_sets, exon_set_cts = exon_set_cts)
  }
  
  return(output)
  
}

#-----------------------------------------------------------------------------#
# Create Pure CT Reference Count Files in Dirichlet-Negative Binomial Model   #
#-----------------------------------------------------------------------------#

counts_output = function(exonInfo_1000, exonInfo_other, theta, file_labels,
                         folder, seed){
  
  set.seed(seed)
  
  # Checks
  
  
  # Define variables
  ## Number cell types
  J = length(exonInfo_1000)
  ## Number samples per cell type
  n = length(theta) / J
  
  output_1000 = list()
  
  for(ct in 1:J){
    
    ct_info = exonInfo_1000[[ct]]
    ct_theta = theta[(1 + n*(ct-1)):(n*ct)]
    
    for(clust in names(ct_info)){
      
      mu = ct_info[[clust]]$mu
      exon_sets = ct_info[[clust]]$exon_sets
      
      # Initialize ES_labels (exon set labels) and record of counts (counts_record) in first cluster
      if(names(ct_info)[1] == clust){
        
        ES_labels = exon_sets
        counts_record = matrix(NA, nrow = length(exon_sets), ncol = n)
        
        for(i in 1:n){
          counts_record[,i] = rnegbin(n = length(mu[,i]), mu = mu[,i], theta = ct_theta[i])
        }
        
      }else{ # End IF
        
        ES_labels = c(ES_labels, exon_sets)
        counts_subset = matrix(NA, nrow = length(exon_sets), ncol = n)
        
        for(i in 1:n){
          counts_subset[,i] = rnegbin(n = length(mu[,i]), mu = mu[,i], theta = ct_theta[i])
        }
        
        counts_record = rbind(counts_record, counts_subset)
        
      } # End ELSE of IF-ELSE
      
    } # End clust for-loop
    
    rownames(counts_record) = ES_labels
    
    output_1000[[ct]] = list(ES_labels = ES_labels, counts = counts_record)
    
  } # End ct for-loop
  
  output_other = list()
  
  for(ct in 1:J){
    
    ct_info = exonInfo_other[[ct]]
    
    for(clust in names(ct_info)){
      
      exon_sets = ct_info[[clust]]$exon_sets
      
      # Initialize ES_labels (exon set labels) and record of counts (counts_record) in first cluster
      if(names(ct_info)[1] == clust){
        
        ES_labels = exon_sets
        counts_record = ct_info[[clust]]$exon_set_cts
        
      }else{ # End IF
        
        ES_labels = c(ES_labels, exon_sets)
        counts_record = rbind(counts_record, ct_info[[clust]]$exon_set_cts)
        
      } # End ELSE of IF-ELSE
      
    } # End clust for-loop
    
    output_other[[ct]] = list(ES_labels = ES_labels, counts = counts_record)
    
  } # End ct for-loop
  
  for(ct in 1:J){
    
    ct_files = file_labels[(1 + n*(ct-1)):(n*ct)]
    
    counts_combo = rbind(output_1000[[ct]]$counts, output_other[[ct]]$counts)
    ES_labels_all = c(output_1000[[ct]]$ES_labels, output_other[[ct]]$ES_labels)
    
    for(i in 1:n){
      df = data.frame(counts = counts_combo[,i], exons = ES_labels_all)
      write.table(df, file = sprintf("%s/%s.txt", folder, ct_files[i]), 
                  row.names = F, col.names = F)
    } # End i for-loop
    
  } # End ct for-loop
  
} # End counts_output

#-----------------------------------------------------------------------------#
# Create Pure CT Reference Count Files in Dirichlet-Multinomial Model         #
#-----------------------------------------------------------------------------#

counts_output_MM = function(exonInfo_1000, exonInfo_other, n, file_labels, folder, seed){
  
  set.seed(seed)
  
  # Checks
  
  # Define variables
  ## Number cell types
  J = length(exonInfo_1000)
  
  output_1000 = list()
  
  for(ct in 1:J){
    
    ct_info = exonInfo_1000[[ct]]
    
    for(clust in names(ct_info)){
      
      gene_ct = ct_info[[clust]]$gene_ct
      mm_probs = ct_info[[clust]]$mm_probs
      exon_sets = ct_info[[clust]]$exon_sets
      
      # Initialize ES_labels (exon set labels) and record of counts (counts_record) in first cluster
      if(names(ct_info)[1] == clust){
        
        ES_labels = exon_sets
        counts_record = matrix(NA, nrow = length(exon_sets), ncol = n)
        
        for(i in 1:n){
          counts_record[,i] = rmultinom(n = 1, size = gene_ct[i], prob = mm_probs[,i])
        }
        
      }else{ # clust not first cluster
        
        ES_labels = c(ES_labels, exon_sets)
        counts_subset = matrix(NA, nrow = length(exon_sets), ncol = n)
        
        for(i in 1:n){
          counts_subset[,i] = rmultinom(n = 1, size = gene_ct[i], prob = mm_probs[,i])
        }
        
        counts_record = rbind(counts_record, counts_subset)
        
      } # End ELSE of IF-ELSE
      
    } # End clust for-loop
    
    rownames(counts_record) = ES_labels
    
    output_1000[[ct]] = list(ES_labels = ES_labels, counts = counts_record)
    
  } # End ct for-loop
  
  output_other = list()
  
  for(ct in 1:J){
    
    ct_info = exonInfo_other[[ct]]
    
    for(clust in names(ct_info)){
      
      exon_sets = ct_info[[clust]]$exon_sets
      
      # Initialize ES_labels (exon set labels) and record of counts (counts_record) in first cluster
      if(names(ct_info)[1] == clust){
        
        ES_labels = exon_sets
        counts_record = ct_info[[clust]]$exon_set_cts
        
      }else{ # End IF
        
        ES_labels = c(ES_labels, exon_sets)
        counts_record = rbind(counts_record, ct_info[[clust]]$exon_set_cts)
        
      } # End ELSE of IF-ELSE
      
    } # End clust for-loop
    
    output_other[[ct]] = list(ES_labels = ES_labels, counts = counts_record)
    
  } # End ct for-loop
  
  for(ct in 1:J){
    
    ct_files = file_labels[(1 + n*(ct-1)):(n*ct)]
    
    counts_combo = rbind(output_1000[[ct]]$counts, output_other[[ct]]$counts)
    ES_labels_all = c(output_1000[[ct]]$ES_labels, output_other[[ct]]$ES_labels)
    
    for(i in 1:n){
      df = data.frame(counts = counts_combo[,i], exons = ES_labels_all)
      write.table(df, file = sprintf("%s/%s.txt", folder, ct_files[i]), 
                  row.names = F, col.names = F)
    } # End i for-loop
    
  } # End ct for-loop
  
} # End counts_output



#-----------------------------------------------------------------------------#
# Simulate Mixture Count Files                                                #
#-----------------------------------------------------------------------------#

mix_creation = function(set_mixSim, out_folder, file_labels, total_cts, probs, seed){
  
  set.seed(seed)
  
  # Define variables
  ## Number mixture replicates to create
  mix_rep = nrow(total_cts)
  ## Number cell types
  J = ncol(probs)
  ## Number pure reference samples per cell type (assume equal across all cell types)
  M = length(set_mixSim[[1]])
  
  # Checks
  if(any(rowSums(probs) != 1)){
    stop("rows of probs must add to 1")
  }
  
  # List of pure reference sample count data.frames
  df_list = list()
  for(j in 1:J){
    pure_files = set_mixSim[[j]]
    files_list = list()
    for(f in 1:length(pure_files)){
      df = read.table(file = pure_files[f], as.is = T)
      colnames(df) = c("count","exons")
      files_list[[f]] = df
    }
    df_list[[j]] = files_list
  }
  
  # exon set labels (assume same across all pure reference samples)
  exon_sets = df_list[[1]][[1]]$exons
  # Number exon sets (assume equal across all pure reference samples)
  E = length(exon_sets)
  
  for(k in 1:nrow(probs)){
    # Identify prob vector
    p = probs[k,]
    # Randomly select counts files from each cell type
    pure_counts = matrix(NA, nrow = E, ncol = J)
    for(j in 1:J){
      counts_vec = df_list[[j]][[sample(1:M, size = 1)]]$count
      pure_counts[,j] = counts_vec
    }
    
    # Calculate ratio of total counts between mixture replicate and pure reference counts
    cts_Ratio = matrix(NA, nrow = mix_rep, ncol = J)
    for(rep in 1:mix_rep){
      cts_Ratio[rep,] = total_cts[rep,k] / colSums(pure_counts)
    }
    
    # Multiply p and cts_Ratio to appropriate columns of pure_counts to get mixture sample components
    ## Round results and add results across exon sets
    mixture = list()
    for(rep in 1:mix_rep){
      mix_components = pure_counts * matrix(p, nrow = nrow(pure_counts), ncol = J, byrow = T) *
        matrix(cts_Ratio[rep,], nrow = nrow(pure_counts), ncol = J, byrow = T)
      mixture[[rep]] = rowSums(round(mix_components))
    }
    
    # Save mixture results in counts.txt files
    for(rep in 1:mix_rep){
      label = file_labels[rep,k]
      df_mix = data.frame(count = mixture[[rep]], exons = exon_sets)
      write.table(df_mix, file = sprintf("%s/%s.txt", out_folder, label), col.names = F, row.names = F)
    }
  }
  
}

mix_creation2 = function(set_mixSim, out_folder, file_labels, total_cts, probs, seed){
  
  set.seed(seed)
  
  # Define variables
  
  ## Number cell types
  J = ncol(probs)
  ## Number pure reference samples per cell type (assume equal across all cell types)
  M = length(set_mixSim[[1]])
  
  # Checks
  if(!all.equal(sum(probs), 1, check.attributes = F)){
    stop("probs must add to 1")
  }
  
  # List of pure reference sample count data.frames
  df_list = list()
  for(j in 1:J){
    pure_files = set_mixSim[[j]]
    files_list = list()
    for(f in 1:length(pure_files)){
      df = read.table(file = pure_files[f], as.is = T)
      colnames(df) = c("count","exons")
      files_list[[f]] = df
    }
    df_list[[j]] = files_list
  }
  
  # exon set labels (assume same across all pure reference samples)
  exon_sets = df_list[[1]][[1]]$exons
  # Number exon sets (assume equal across all pure reference samples)
  E = length(exon_sets)
  
  for(k in 1:nrow(probs)){
    # Identify prob vector
    p = probs[k,]
    # Randomly select counts files from each cell type
    pure_counts = matrix(NA, nrow = E, ncol = J)
    for(j in 1:J){
      counts_vec = df_list[[j]][[sample(1:M, size = 1)]]$count
      pure_counts[,j] = counts_vec
    }

    # Calculate ratio of total counts between mixture sample and pure reference counts
    # Goal: standardize total counts from each pure sample, then take desired proportion
    cts_Ratio = total_cts[k] / colSums(pure_counts)
    
    # Multiply p and cts_Ratio to appropriate columns of pure_counts to get mixture sample components
    mix_components = pure_counts * matrix(p, nrow = nrow(pure_counts), ncol = J, byrow = T) *
      matrix(cts_Ratio, nrow = nrow(pure_counts), ncol = J, byrow = T)
    mixture = rowSums(round(mix_components))

   # Save mixture results in counts.txt files
    label = file_labels[k]
    df_mix = data.frame(count = mixture, exons = exon_sets)
    write.table(df_mix, file = sprintf("%s/%s.txt", out_folder, label), col.names = F, row.names = F)
    
  }
  
  # # Identify prob vector
  # p = probs
  # # Randomly select counts files from each cell type
  # pure_counts = matrix(NA, nrow = E, ncol = J)
  # for(j in 1:J){
  #   counts_vec = df_list[[j]][[sample(1:M, size = 1)]]$count
  #   pure_counts[,j] = counts_vec
  # }
  # 
  # # Calculate ratio of total counts between mixture replicate and pure reference counts
  # cts_Ratio = matrix(NA, nrow = mix_rep, ncol = J)
  # for(rep in 1:mix_rep){
  #   cts_Ratio[rep,] = total_cts[rep,k] / colSums(pure_counts)
  # }
  # 
  # # Multiply p and cts_Ratio to appropriate columns of pure_counts to get mixture sample components
  # ## Round results and add results across exon sets
  # mixture = list()
  # for(rep in 1:mix_rep){
  #   mix_components = pure_counts * matrix(p, nrow = nrow(pure_counts), ncol = J, byrow = T) *
  #     matrix(cts_Ratio[rep,], nrow = nrow(pure_counts), ncol = J, byrow = T)
  #   mixture[[rep]] = rowSums(round(mix_components))
  # }
  # 
  # # Save mixture results in counts.txt files
  # for(rep in 1:mix_rep){
  #   label = file_labels[rep,k]
  #   df_mix = data.frame(count = mixture[[rep]], exons = exon_sets)
  #   write.table(df_mix, file = sprintf("%s/%s.txt", out_folder, label), col.names = F, row.names = F)
  # }

  
}

#-----------------------------------------------------------------------------#
# Simulate Fragment Length Distribution Files                                 #
#-----------------------------------------------------------------------------#

fragLens_out = function(total_reads, mean = 300, SD = 50, lenMin = 150, lenMax = 600,
                        out_folder, file_names, seed){
  
  set.seed(seed)
  
  # Define variables
  mix_rep = nrow(total_reads)
  num_pCombos = ncol(total_reads)
  
  for(rep in 1:mix_rep){
    
    for(p in 1:num_pCombos){
      # fragLens_dist() in geneModel_code/fragLens_dist.cpp file
      freq_dist = fragLens_dist(total_reads[rep,p], mean, SD, lenMin, lenMax)
      freq_dist = freq_dist[which(freq_dist[,1] > 0),]
      write.table(freq_dist, file = sprintf("%s/%s.txt",out_folder,file_names[rep,p]), col.names = F, row.names = F)
    }
    
  }
  
}
