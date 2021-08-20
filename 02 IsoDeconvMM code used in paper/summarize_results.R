#-----------------------------------------------------------------------#
# SUMMARIZING the Output                                                #
#-----------------------------------------------------------------------#
# Compile the estimates from each cluster, provide a histogram, 
# and summarize the values with a geometric median.

#' @importFrom ICSNP spatial.median
#' @import ggplot2
#' @export
Summarize_Report<-function(final_output, restrict = FALSE, clusters = NULL, 
                           plots_options = plotsConrol()){
  
  mix_names = names(final_output)
  
  if(restrict==TRUE && length(clusters)>0){
    warn("Set restrict to FALSE, but specified clusters for restriction.\n")
    warn("Will restrict to provided clusters.\n")
  }
  
  final_summary = list()
  
  for(j in 1:length(final_output)){
    cdata = final_output[[j]]
    
    #--------------------------------------------------------------#
    # CODING for restrictions                                      #
    #--------------------------------------------------------------#
    if(length(clusters)>0){
      cdata_fin = cdata[clusters]
    } else {
      cdata_fin = cdata
    }
    
    #--------------------------------------------------------------#
    # EXTRACT cell-types info                                      #
    #--------------------------------------------------------------#
    cto = cdata_fin[[1]][["CellType_Order"]]
    nclust = length(cdata_fin)
    clust_names = names(cdata_fin)
    
    p_mat = matrix(0, nrow = nclust, ncol = length(cto))
    colnames(p_mat) = cto
    rownames(p_mat) = clust_names
    warn = numeric(nclust)
    
    for(i in 1:nclust){
      p_mat[i,] = cdata_fin[[i]][["mix"]][["p.est"]]
      warn[i] = cdata_fin[[i]][["WARN"]]
    }
    
    if(any(p_mat< -1e-5)){message("Warning: Prop < -1e-5")}
    if(any(p_mat>(1+1e-5))){message("Warning: Prop > 1 + 1e-5")}
    
    p_mat[which(p_mat<0)] = 0
    # p_mat[which(p_mat>1)] = 1
    
    #--------------------------------------------------------------#
    # PLOTTING the cell type info                                  #
    #--------------------------------------------------------------#
    
    if(plots_options$plots == TRUE){
      df = as.data.frame(p_mat)
      
      histograms = lapply(1:2, function(k){
        mtitle = sprintf("Distribution of Estimated Proportions of (%s)",cto[k])
        xlabel= sprintf("Prop. of (%s)",cto[k])
        hg = ggplot(data = df, mapping = aes(df[,k])) + geom_histogram(bins = plots_options$bins) +
          ggtitle(label = mtitle, subtitle = sprintf("Mixture %s", mix_names[j])) + xlab(xlabel)
      })
      
      names(histograms) = cto
    }
    
    #---------------------------------------------------------------#
    # Summarizing the Values                                        #
    #---------------------------------------------------------------#
    fin_est = matrix(0,nrow=1,ncol=length(cto))
    colnames(fin_est) = cto
    
    #------------------------------------------------------#
    # WARNING INDICATORS                                   #
    #------------------------------------------------------#
    # 0 - Optimization Complete
    # 1 - Iteration Limit Reached
    # 4 - Error in Optimization Routine
    # 5 - Optimization not conducted (Error in pure sample fit)
    # Note: Should only evaluate p estimates where WARN is only 0 or 1 (?)
    # p.fin = spatial.median(X = p_mat[,-length(cto)])
    # p.fin = spatial.median(X = p_mat[which(warn <= 1),-length(cto)])
    # 
    # fin_est[1,] = c(p.fin,1-sum(p.fin))
    
    # if(plots_options$plots == TRUE){
    #   final_summary[[mix_names[j]]] = list(p_est = fin_est, p_mat = data.frame(p_mat, WARN = warn), histograms = histograms)
    # }else{
    #   final_summary[[mix_names[j]]] = list(p_est = fin_est, p_mat = data.frame(p_mat, WARN = warn))
    # }
    
    if(plots_options$plots == TRUE){
      final_summary[[mix_names[j]]] = list(p_mat = data.frame(p_mat, WARN = warn), histograms = histograms)
    }else{
      final_summary[[mix_names[j]]] = list(p_mat = data.frame(p_mat, WARN = warn))
    }
    
  }
  
  #----------------------------------------------------------------#
  # RETURN values                                                  #
  #----------------------------------------------------------------#
  return(final_summary)
  
}

#' @export
plotsControl = function(plots = TRUE, bins = 20){
  structure(list(plots = plots, bins = bins), 
            class = c("plotsControl"))
}