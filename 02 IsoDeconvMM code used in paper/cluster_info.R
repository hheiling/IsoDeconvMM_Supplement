
#' @importFrom stringr str_split str_detect
cluster_info = function(countData, labels, cellTypes, bedFile, discrim_genes, discrim_clusts,
                        readLen, lmax){
  
  #------------------------------------------------------------------------------------------------------------------------------------#
  # LOADING THE DATA                                                                                                                   #
  #------------------------------------------------------------------------------------------------------------------------------------#
  
  # Edit countData to only include rows with the discrim_genes or discrim_clusts
  row_remove = function(counts, discrim_genes, discrim_clusts){
    # Restrict rows of count data data.frame based on discrim_genes or discrim_clusts
    if(!is.null(discrim_clusts)){
      discrim_units = discrim_clusts
    }else if(!is.null(discrim_genes)){
      discrim_units = discrim_genes
    }else if(is.null(discrim_genes) & is.null(discrim_clusts)){
      return(counts)
    }
    
    row_idx = c()
    for(u in 1:length(discrim_units)){
      row_idx = c(row_idx, which(str_detect(as.character(counts$exons),as.character(discrim_units[u]))))
    }
    
    counts_new = counts[row_idx,]
    return(counts_new)
  }
  
  # Calls to Dr. Sun's loadData function and generates one GeneModel for each sample
  # present in the countData vector.
  
  # See "R/loadData_edit.R" file for loadData_djEdit() function code
  for(i in 1:length(countData)){
    # ct_datai = countData[[i]]
    # Restrict ct_datai to only contain data from clusters of interest?
    ## Will this have any negative downstream effects?
    ct_datai = row_remove(countData[[i]], discrim_genes, discrim_clusts)
    outFile = labels[i]
    cat("\n Loading data for ", outFile, "\n")
    cmdi = sprintf("%s = loadData_djEdit(ct_datai,bedFile,readLen,lmax)",outFile)
    eval(parse(text=cmdi))
  }
  
  #------------------------------------------------------------------------------------------------------------------------------------#
  # LIST APPEND: ALL SAMPLES                                                                                                           #
  #------------------------------------------------------------------------------------------------------------------------------------#
  
  sam_names = labels
  
  # list_append function: concatenates all geneModels into a large list, labeled by sample names.
  
  list_append <- function(list_names){
    nl = length(list_names)
    concat_list = list()
    for(i in 1:nl){
      cmdi = sprintf("concat_list[[i]] = %s",list_names[i])
      eval(parse(text=cmdi))
    }
    names(concat_list) = list_names
    return(concat_list)
  }
  
  concat_list = list_append(sam_names)
  
  #------------------------------------------------------------------------------------------------------------------------------------#
  # OBTAINING LIST OF TRANSCRIPT CLUSTER NAMES                                                                                         #
  #------------------------------------------------------------------------------------------------------------------------------------#
  
  tnames = NULL
  
  for(i in 1:length(sam_names)){
    call1 = sprintf("tclust_i = names(%s)",sam_names[i])
    eval(parse(text=call1))
    call2 = sprintf("tclust_iand1 = names(%s)",sam_names[i+1])
    eval(parse(text=call2))
    if(i==1){tnames = union(tclust_i,tclust_iand1)}
    else{tnames = union(tnames,tclust_i)}
  }
  
  tnames = tnames[order(tnames)]
  
  # tnames_discrim = tnames[which(tnames %in% discrim_clusters)]
  
  # Create a list of length length(tnames) wherein the i-th element of the list is a list composed of #samples+1 datasets and a row vector.
  # The vector is info_status and confirms whether or not the info dataset is complete and accurate.
  # The first dataset is $info, which contains information on all exons in the cluster. The subsequent datasets are the counts at each
  # exon set in each sample.
  
  concat_geneMod = list()
  
  for(i in 1:length(tnames)){
    concat_geneMod[[i]] = list(info_status="Not Checked",info=NULL)
    incurr_clust = rep(0,length(sam_names))
    z = rep(TRUE,length(sam_names)-1)
    for(j in 1:length(sam_names)){
      tclust_currsamp = names(concat_list[[j]])
      if(tnames[i] %in% tclust_currsamp){
        cmdij = sprintf("concat_geneMod[[i]]$%s = concat_list[[j]][[tnames[i]]]$count",sam_names[j])
        eval(parse(text=cmdij))
        incurr_clust[j] = 1
      } else {
        cmdij = sprintf("concat_geneMod[[i]]$%s = data.frame(count = as.numeric(NULL),exons=as.character(NULL))",sam_names[j])
        eval(parse(text=cmdij))
        incurr_clust[j] = 0
      }
    }
    sam_clust = which(incurr_clust==1)
    fsamp = sam_clust[1]
    concat_geneMod[[i]]$info = concat_list[[fsamp]][[tnames[i]]]$info
    for(k in 1:length(sam_clust)){
      fsamp_k = sam_clust[k]
      dim_orig = dim(concat_geneMod[[i]]$info)
      dim_final = dim(concat_list[[fsamp_k]][[tnames[i]]]$info)
      if(all(dim_orig == dim_final)){
        z[k-1] = all(concat_geneMod[[i]]$info==concat_list[[fsamp_k]][[tnames[i]]]$info)
      } else {
        z[k-1] = FALSE
      }
    }
    if (all(z)){
      concat_geneMod[[i]]$info_status = "PASS"
    } else {
      concat_geneMod[[i]]$info_status = "FAIL"
    }
  }
  
  names(concat_geneMod) = tnames
  
  #-------------------------------------------------------------------------------------------------------------------------------------------#
  # CONCAT_GENEMOD:                                                                                                                           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # FORMAT:                                                                                                                                   #
  #     List                                                                                                                                  #
  # I-TH ELEMENT:                                                                                                                             #
  #     - Label     : "Transcript Cluster I"                                                                                                  #                                                                                                   #
  #     - Element 1 : $info_status (vector)                                                                                                   #
  #         - "Not Checked" -- error in the program and $info data set was not checked for correctness.                                       #
  #         - "PASS" -- $info dataset passes checks for accuracy                                                                              #
  #         - "FAIL" -- $info dataset failed checks for accuracy                                                                              # 
  #     - Element 2 : $info dataset                                                                                                           #
  #         - chr   : contains the chromosome location of the exon being considered                                                           #
  #         - start : contains the starting location of the exon being considered                                                             #
  #         - end   : contains the ending location of the exon being considered                                                               #
  #         - exon  : contains the unique exon id for the exon being considered                                                               #
  #     - Element 3 : $sample_name dataset                                                                                                    #   
  #         - exons : Contains the IDs of the exons comprising the exon set                                                                   #
  #         - count : contains the number of reads at the exon set in question for "sample_name"                                              #
  #-------------------------------------------------------------------------------------------------------------------------------------------#
  
  #-----------------------------------------------------------------------------#
  # ESTABLISHING CLUSTERS WITH HIGHEST LEVELS OF DISCRIMINATORY CAPABILITY      #
  #-----------------------------------------------------------------------------#
  # Limit clusters examined to those with discriminatory genes 
  # (discrim_genes or discrim_clusts)                                                                                                                       #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Gene names for cluster given in "info" matrix                               #
  #-----------------------------------------------------------------------------#
  
  #------------------ Identify Highly Discriminatory Clusters -------------------#
  
  # User inputs discrim_genes or discrim_clusts information
  # Find names of clusters that contain these discriminatory genes
  
  if(!is.null(discrim_clusts)){
    
    discrim_clusters = discrim_clusts
    
  }else if(!is.null(discrim_genes)){
    
    all_clusters = names(concat_geneMod)
    idx_clust_tmp = numeric(length(all_clusters))
    
    
    for(clust in all_clusters){
      clust_genes = unique(unlist(str_split(concat_geneMod[[clust]]$info$gene, pattern = ":")))
      if(any(clust_genes %in% discrim_genes)){
        idx_clust_tmp[which(all_clusters == clust)] = 1
      }
    }
    
    idx_clust = which(idx_clust_tmp==1)
    discrim_clusters = unique(all_clusters[idx_clust])
    
  }
  
  cat("Length discrim_clusters: ", length(discrim_clusters), "\n")
  
  concat_geneMod = concat_geneMod[discrim_clusters]
  
  return(concat_geneMod)
  
}

