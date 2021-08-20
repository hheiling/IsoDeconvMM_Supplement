

#' @title Cell Type Deconvolution using RNA Isoform-Level Expression
#'  
#' @description Calculates the proportions of pure cell type components in heterogeneous cell type samples of 
#' RNA-seq data utilizing isoform-level expression differences
#' 
#' @param directory an optional character string denoting the path to the directory where all of the 
#' mix_files, pure_ref_files, fraglens_files, and bedfile are located. The working directory is set as 
#' this directory. If this directory is left `NULL`, then all of the relevent files must either (a) be 
#' located in the current working directory or (b) have their full path specified.
#' @param mix_files a vector of the file names for the text files recording the number of RNA-seq 
#' fragments per exon set, which should  have 2 columns "count" and "exons", without header. 
#' For example:
#' \verb{
#'    37  chr18_109|ENSMUSG00000024491|4;chr18_109|ENSMUSG00000024491|5;
#'    17  chr18_109|ENSMUSG00000024491|5;
#'    88  chr18_109|ENSMUSG00000024491|5;chr18_109|ENSMUSG00000024491|6;
#' }
#' There should be one file for each of the samples containing mixtures of cells. 
#' The second column lists exon sets, where “chr18_109” indicates a transcript cluster, 
#' “ENSMUSG00000024491” is the ensemble gene ID, and the numbers at the end is the exon ID
#' Directions to create these count files can be found in the Step_0_Processes directory of the GitHub repo 
#' hheiling/deconvolution <https://github.com/hheiling/deconvolution>
#' @param pure_ref_files a matrix where the first column is the file names for the text files
#' recording the number of RNA-seq fragments per exon set (see `mix_files` for additional description),
#' one for each of the pure reference cell type samples (again, see the Step_0_Processes directory in
#' <https://github.com/hheiling/deconvolution> for directions on how to create these files) and the 
#' second column contains the character names of the pure cell type associated with each sample
#' @param fraglens_files a vector of the file names for the text files recording the distribution 
#' of the fragment lengths, which should have 2 columns: "Frequency" and "Length", without header.
#' For example:
#' \verb{
#'    20546 75
#'    40465 76
#'    37486 77
#'    27533 78
#'    25344 79
#' }
#' Directions to create these fragment length files are also available in the Step_0_Processes directory in
#' the GitHub repo hheiling/deconvoltuion, <https://github.com/hheiling/deconvolution>
#' @param bedFile file name of the .bed file recording information of non-overlapping exons, which  
#' has 6 colums: "chr", "start", "end", "exon", "score", and "strand",
#' without header. For example:
#' \verb{
#'   chr1    3044314 3044814 ENSMUSG00000090025:1    666     +
#'   chr1    3092097 3092206 ENSMUSG00000064842:1    666     +
#' }
#' Directions to create this .bed file can be found in the 
#' Create_BED_knownIsoforms_Files directory in the GitHub repo hheiling/deconvolution, 
#' <https://github.com/hheiling/deconvolution>
#' @param knownIsoforms character string for the name of an .RData object that contains the known isoform
#' information. When loaded, this object is a list where each component is a binary matrix 
#' that specifies a set of possible isoforms (e.g., isoforms from annotations). Specifically, it is a 
#' binary matrix of k rows and m columns, where k is the number of 
#' non-overlapping exons and m is the number of isoforms. isoforms[i,j]=1 
#' indicates that the i-th exon belongs to the j-th isoform. For example, 
#' the following matrix indicates the three isoforms for one gene ENSMUSG00000000003:
#' \verb{
#'      ENSMUST00000000003 ENSMUST00000166366 ENSMUST00000114041
#' [1,]                  1                  1                  1
#' [2,]                  1                  1                  1
#' [3,]                  1                  1                  1
#' [4,]                  1                  1                  0
#' [5,]                  1                  1                  1
#' [6,]                  1                  1                  1
#' [7,]                  1                  1                  1
#' [8,]                  1                  0                  0
#' }
#' Instructions for creating such an RData object can be found in the
#' Create_BED_knownIsoforms_Files directory in the GitHub repo hheiling/deconvolution, 
#' <https://github.com/hheiling/deconvolution>
#' @param discrim_genes vector of genes that are suspected to have differential gene expression. 
#' This gene list could come from CuffLinks output, \code{isoform} package output, or something similar.
#' @param readLen numeric value of the length of a read in the RNAseq experiment
#' @param lmax numeric value of the maximum fragment length of the experiment
#' @param eLenMin numeric value of the minimum value of effective length. 
#' If the effective length of an exon or exon junction is smaller than eLenMin,
#' i.e., if this exon is not included in the corresponding isoform, 
#' set it to eLenMin. This is to account for possible sequencing error or
#' mapping errors.
#' @param mix_names an optional vector of the desired nicknames of the mixture samples corresponding,
#' in the same order, to the mix_files list. If left as the default \code{NULL} value, the nicknames used
#' will be the names given in the mix_files minus the .txt extension
#' @param initPts an optional matrix of initial probability estimates for the cell composition
#' of the mixture samples to be used in the optimization procedure. The matrix should have J columns,
#' where J = number of pure cell types of interest. Each row corresponds to different combinations of 
#' initial probability values. The column names of the matrix must be provided and must correspond 
#' to the pure cell type names given in the second column of the pure_ref_files object (no particular ordering needed)
#' @param optim_options a list inheriting from class \code{optimControl} containing optimization 
#' control parameters. See the function \code{\link{optimControl}} for more details.
#' 
#' @return A list object with the following structure: first layer of list has elements associated with each
#' of the mixture samples; second layer of list as elements associated with each transcript cluster
#' used in the analysis, determined by the genes in the discrim_gene vector. Each of these transcript 
#' cluster elements is itself a list with the following elements:
#' \item{info}{}
#' \item{candiIsoform}{}
#' \item{I}{Number of isoforms utilized in transcript cluster}
#' \item{E}{Number of exons in transcript cluster}
#' \item{X}{ExI matrix of effective lengths for each of the E exon sets within each of the I isoforms}
#' \item{info_status}{}
#' \item{y_mix, other y vectors for each pure cell type reference sample}{Ex1 vectors of read count at 
#' each exon set for the given mixture or pure cell type sample}
#' \item{countN_mix, other countN values for each pure cell type reference sample}{}
#' \item{mix}{a list with the elements rds_exons_t (vector of length E+1 where the last E elements are y_mix,
#' and the first element is the total read counts for the mixture sample minus the sum of y_mix), 
#' gamma.est ((I-1)xK matrix of isoform expression parameters for each cell type k), 
#' tau.est (vector of length K of gene expression parameters in cell type k),
#' p.est (vector of length K containing estimated proportions based on the given transcript cluster),
#' and pm.rds.exons (ExK matrix containing posterior means for each of E exon sets in each of K cell types)
#' }
#' \item{"cellType1","cellType2" ...}{}
#' \item{l_tilde}{Ix1 vector of total effective lengths of each of the I isoforms; 
#' Each elemement of the vector, denoted l_i, is a column sum from the matrix \code{X}}
#' \item{X.fin}{edited design matrix for new gamma parameters, where the ith column of the new matrix is
#' \verb{X.fin_i = (X_i-[l_i/l_I]X_I)} for i = 1,...,(I-1) and \verb{X.fin_I = X_I/l_I}}
#' \item{X.prime}{first (I-1) columns X.fin pertaining to gamma parameters}
#' \item{alpha.est}{IxK hyperparameters governing average isoform expression levels and variances within cells of type k}
#' \item{beta.est}{2xK hyperparameters governing gene expression levels within cells type k}
#' \item{CellType_Order}{For outputs giving K different estimates for each of the K cell types,
#' these outputs are ordered with respect to CellType_Order}
#' \item{WARN}{An integer indicating the following information:
#' \verb{
#'              0 - Optimization Complete
#'              1 - Iteration Limit Reached
#'              4 - Error in Optimization Routine (Error in mixture sample fit)
#'              5 - Optimization not conducted (Error in pure sample fit)
#'              }}
#'  
#' @importFrom stringr str_c str_remove
#' @export
IsoDeconvMM = function(directory = NULL, mix_files, pure_ref_files,
                       fraglens_files, bedFile, knownIsoforms, 
                       discrim_genes = NULL, discrim_clusts = NULL,
                       readLen, lmax = 600, eLenMin = 1, mix_names = NULL,
                       initPts = NULL,
                       optim_options = optimControl()){
  
  time0 = proc.time()
  
  if(is.null(discrim_genes) & is.null(discrim_clusts)){
    stop("Either discrim_genes or discrim_clusts must be specified \n")
  }else if(!is.null(discrim_genes) & !is.null(discrim_clusts)){
    warning("Both discrim_genes and discrim_clusts specified, ",
            "will only use discrim_clusts \n", immediate. = T)
  }
  
  #-----------------------------------------------------------------------------#
  # Download all needed files, convert input items to useful formats            #
  #-----------------------------------------------------------------------------#
  
  if(!is.null(directory)){
    setwd(directory)
  }
  
  countData_mix = mix_files
  countData_pure = pure_ref_files[,1] 
  
  cellTypes_pure = as.character(pure_ref_files[,2])
  # ctpure_names = pure cellTypes names
  ctpure_names = levels(as.factor(cellTypes_pure))
  
  labels_pure = character(length(countData_pure))
  
  for(type in ctpure_names){
    locations = which(cellTypes_pure == type)
    labels_type = cellTypes_pure[locations]
    num_ref = length(labels_type)
    labels_ref = str_c(labels_type, "_ref", 1:num_ref)
    labels_pure[locations] = labels_ref
  }
  
  if(!is.null(mix_names)){
    labels_mix = mix_names
  }else{
    labels_mix = str_c("mix",1:length(mix_files))
  }
  
  # Download all count text files, compute total counts for each file
  # Output: list with elements total_cts, counts_list
  ## See comp_total_cts() function under "Internal isoDeconvMM Functions" heading later in this document
  pure_input = comp_total_cts(directory = directory, countData = countData_pure)
  pure_counts = pure_input$counts_list
  
  mix_input = comp_total_cts(directory = directory, countData = countData_mix)
  mix_counts = mix_input$counts_list
  
  countData = pure_counts
  countData[(length(pure_counts)+1):(length(pure_counts)+length(mix_counts))] = mix_counts
  
  if(!(length(fraglens_files) %in% c(1, length(mix_files)))){
    stop("Length of fraglens_files must be equal to 1 or length of mix_files")
  }
  
  fraglens_list = list()
  for(i in 1:length(fraglens_files)){
    fraglens = read.table(fraglens_files[i], as.is = T)
    if (ncol(fraglens) != 2) {
      stop(fraglens_files[i], " should have 2 columns for Freq and Len\n")
    }
    names(fraglens) = c("Freq", "Len")
    
    fraglens_list[[i]] = fraglens
  }
  
  
  # CHECKING Presence of Isoforms File:
  # IsoDeconv requires a list of known isoforms in order to model intra-sample heterogeneity. Stops program if file not present.
  # Load knownIsoforms .RData object:
  
  if(!is.null(knownIsoforms)){
    assign("isoAll", get(load(knownIsoforms)))
  }else{
    stop("knownIsoforms list object must be present!")
  }
  
  # Download .bed file
  bedFile_info = read.table(sprintf("%s", bedFile), sep = "\t", as.is = TRUE)
  
  bf_colNames = c("chr", "start", "end", "exon", "score", "strand")
  
  if (ncol(bedFile_info) != 6) {
    cN = paste(bf_colNames, collapse = ", ")
    stop(bedFile, " should have 6 columns: ", cN, "\n")
  }
  
  names(bedFile_info) = bf_colNames
  
  print("Finished loading supporting files")
  time1 = proc.time()
  
  # Check initPts is specified correctly OR create initPts matrix if not specified
  
  if(is.null(initPts)){
    ct_num = length(ctpure_names)
    if(ct_num == 2){
      p1 = c(0.1, 0.25, 0.33, 0.5, 0.67, 0.75, 0.9)
      initPts = cbind(p1, 1-p1)
    }else if(ct_num == 3){
      initPts = matrix(c(0.10,0.10,0.80,
                         0.10,0.80,0.10,
                         0.80,0.10,0.10,
                         0.25,0.25,0.50,
                         0.25,0.50,0.25,
                         0.50,0.25,0.25,
                         0.20,0.40,0.40,
                         0.40,0.20,0.40,
                         0.40,0.40,0.20,
                         1/3,1/3,1/3), ncol = 3, byrow = T)
    }else if((ct_num >= 4) & (ct_num <= 9)){
      initPts = matrix(0.1, nrow = ct_num, ncol = ct_num)
      for(j in 1:ct_num){
        initPts[j,j] = 1 - 0.1 * ct_num
      }
      initPts = rbind(initPts, rep(1/ct_num, times = ct_num))
    }else if(ct_num >= 10){
      initPts = matrix(1/ct_num, nrow = 1, ncol = ct_num)
    }
    
    # initPts = matrix(1/length(ctpure_names), nrow = 1, ncol = length(ctpure_names))
    colnames(initPts) = ctpure_names
  }else if(class(initPts) == "matrix"){
    if(ncol(initPts) != length(ctpure_names)){
      stop("number of columns of initPts must be equal to the number of pure reference cell types")
    }
    if(!all(colnames(initPts) %in% ctpure_names)){
      stop("colnames of initPts must match pure cell type names in pure_ref_files; see documentation for more details")
    }
    if(sum(initPts > 1) > 0 | sum(initPts < 0) > 0 | (nrow(initPts) > 1 & sum(rowSums(initPts) > 1) > 0)){
      stop("initPts specifies probabilities, which cannot be below 0 or above 1 or sum to more than 1")
    }
    
    # Re-arrange column order of initPts matrix to match order of ctpure_names
    initPts = initPts[,ctpure_names]
    if(!is.matrix(initPts)){
      # If only one initial point, convert above result from vector to matrix with one row
      initPts = matrix(initPts, nrow = 1)
    }
    colnames(initPts) = ctpure_names
    
  }else{
    stop("initPts must be a matrix, see documentation for details")
  }
  
  #--------------------------------------------------------------------------------#
  # Step 1
  # concat_geneMod() extracts gene info, counts, and exon information
  # geneModel_creation() calls geneModel_multcell_Edit(), which is an edit of Dr. Wei Sun's
  # geneModel creation problem (after edits, now accommodates multiple cell types)
  # geneModel() altered by:                                                                                                                        
  #    Douglas Roy Wilson, Jr. 
  #--------------------------------------------------------------------------------#
  
  labels = c(labels_pure, labels_mix)
  cellTypes = c(cellTypes_pure, rep("mix", times = length(mix_files)))
  
  concat_geneMod = cluster_info(countData = countData, labels = labels, cellTypes = cellTypes, 
                                bedFile = bedFile_info, discrim_genes = discrim_genes,
                                discrim_clusts = discrim_clusts, readLen = readLen, lmax = lmax)
  
  final_geneMod = list()
  
  for(j in 1:length(mix_files)){
    cellTypes_sub = c(cellTypes_pure, "mix")
    labels_sub = c(labels_pure, labels_mix[j])
    total_cts = c(pure_input$total_cts, mix_input$total_cts[j])
    if(length(fraglens_list) == 1){
      fragSizeFile = fraglens_list[[1]]
    }else if(length(fraglens_list) == length(mix_files)){
      fragSizeFile = fraglens_list[[j]]
    }
    
    # Need to select specific components of concat_geneMod?
    
    fin_geneMod = geneModel_creation(concat_geneMod = concat_geneMod, fragSizeFile = fragSizeFile, 
                                     cellTypes = cellTypes_sub, labels = labels_sub,
                                     knownIsoforms = isoAll, readLen = readLen, lmax = lmax, 
                                     eLenMin = eLenMin, total_cts = total_cts)
    
    # Perform some checks on the geneMod output
    ## See R/rem_clust.R for rem_clust() code
    sig_geneMod = rem_clust(geneMod = fin_geneMod,co = 5,min_ind = 0)
    
    final_geneMod[[j]] = sig_geneMod
    
  }
  
  print("Finished creation of gene model")
  time2 = proc.time()
  
  #--------------------------------------------------------------------------------#
  # Step 2
  # Add rds_exons matrix for each cell type to each cluster element
  # rds_exons matrix: columns = samples associated with each cell type,
  # rows (except first) = read count for each exon set in the given gene/cluster 
  # for sample j of cell type k,
  # first row = total read count outside gene/cluster of interest in sample j of 
  # cell type k
  # EDIT TO GROUP CELL TYPES
  #--------------------------------------------------------------------------------#
  
  ## See mod_sig_gM() function under "Internal isoDeconvMM Functions" heading later in this document
  modified_sig_geneMod = mod_sig_gM(significant_geneMod = final_geneMod)
  
  #-----------------------------------------------------------#
  # Step 3
  # Pure Cell Type Parameter Estimation                       
  #-----------------------------------------------------------#
  
  ## See pure_estimation() function under "Internal isoDeconvMM Functions" heading later in this document
  pure_est = pure_estimation(modified_sig_geneMod = modified_sig_geneMod, cellTypes = ctpure_names)
  
  print("Finished pure cell type parameter estimation")
  time3 = proc.time()
  
  # Step 4
  #------------------------------------------------------------------------------#
  # Step 4
  # Mixture Cell Type Parameter Estimation
  # Calls STG.Update_Cluster.All() for this estimation procedure
  #------------------------------------------------------------------------------#
  
  IsoDeconv_Output = list()
  
  for(i in 1:length(pure_est)){
    
    tmp.data = pure_est[[i]]
    
    # Establish input break ups
    
    # Data Set Necessities:
    clust.start = 1
    clust.end = length(tmp.data)
    by.value = 15
    
    start.pts = seq(from = 1,to = clust.end,by = by.value)
    end.pts = c((start.pts[-1]-1),clust.end)
    
    cluster_output = list()
    for(m in 1:length(start.pts)){
      start.pt = start.pts[m]
      end.pt = end.pts[m]
      
      curr.clust.opt = tmp.data[c(start.pt:end.pt)]
      ## See R/Production_Functions_MixedSamp.R for STG.Updat_Cluster.All() code
      curr.clust.out = STG.Update_Cluster.All(all_data=curr.clust.opt, cellTypes = ctpure_names,
                                              optimType=optim_options$optimType, 
                                              simple.Init=optim_options$simple.Init, 
                                              initPts = initPts)
      
      cluster_output[[m]] = curr.clust.out
    }
    
    IsoDeconv_Output[[i]] = cluster_output
    
    cat("Finished mixture param estimation for mix file ", j, "\n")
  }
  
  print("Finished mixture paramter estimation for all samples")
  time4 = proc.time()
  
  #---------------------------------------------------------------------------------------------#
  # Step 5
  # Re-compile Step 4 output such that all information organized as follows:
  # First layer of Final_Compiled_Output list is associated with a mixture file of 
  # Second layer of list is associated with name of a cluster
  # Third layer of list contains all cluster-specific information
  #---------------------------------------------------------------------------------------------#
  
  Final_Compiled_Output = list()
  
  for(j in 1:length(IsoDeconv_Output)){
    
    comp.out= NULL
    curr.clust.out = NULL
    
    #---- Set up new pattern ----#
    est.chunks = IsoDeconv_Output[[j]]
    
    message("File ", j)
    
    #---- Populate Output Dataset ----#
    comp.out = list()
    
    for(i in 1:length(est.chunks)){
      curr.clust.out = est.chunks[[i]]
      clust_names = names(curr.clust.out)
      nl = length(curr.clust.out)
      for(m in 1:nl){
        comp.out[[clust_names[m]]]=curr.clust.out[[m]]
      }
    }
    
    Final_Compiled_Output[[j]] = comp.out 
    
  }
  
  if(is.null(mix_names)){
    names(Final_Compiled_Output) = str_remove(mix_files, ".txt")
  }else{
    names(Final_Compiled_Output) = mix_names
  }
  
  # Output timings
  time_record = rbind(time0[1:3],time1[1:3],time2[1:3],time3[1:3],time4[1:3])
  rownames(time_record) = c("Start","End Loading","End Gene Model","End Pure Fit","End Mixture Fit")

  Final_Compiled_Output[["Time"]] = time_record
  
  return(Final_Compiled_Output)
  
  
} # End isoDeconvMM() function



#-------------------------------------------------------------------#
# Internal isoDeconvMM functions                                    #
#-------------------------------------------------------------------#

#-----------------------------------------------------------------------#
# Step 2
# Add rds_exons matrix for each cell type to each cluster list element 
#-----------------------------------------------------------------------#
mod_sig_gM = function(significant_geneMod){
  
  modified_sig_geneMod = list()
  
  for(f in 1:length(significant_geneMod)){
    
    sig_geneMod = significant_geneMod[[f]]
    
    info_mat = sig_geneMod[["Sample_Info"]]
    cellTypes = unique(info_mat$Cell_Type)
    
    ctList = list()
    
    for(j in 1:length(cellTypes)){
      idx = which(info_mat$Cell_Type==cellTypes[j])
      ctList[[cellTypes[j]]] = list(samps = info_mat$Label[idx], tots = info_mat$Total[idx])
    }
    
    idx2consider = which(names(sig_geneMod)!="Sample_Info")
    for(k in idx2consider){
      for(l in 1:length(cellTypes)){
        samps2use = ctList[[l]]$samps
        tots      = ctList[[l]]$tots
        
        y_vecs  = paste("sig_geneMod[[k]]$y",samps2use,sep = "_")
        y_vecsc = paste(y_vecs,collapse = ",")
        nExon = eval(parse(text=sprintf("length(%s)",y_vecs[1])))
        textcmd = sprintf("matrix(c(%s),nrow=nExon,ncol=length(samps2use))",y_vecsc)
        expMat  = eval(parse(text=textcmd))
        
        totmg   = tots - colSums(expMat)
        expMat2 = rbind(totmg,expMat)
        
        if(cellTypes[l]!="mix"){
          sig_geneMod[[k]][[cellTypes[l]]] = list(cellType=cellTypes[l],rds_exons=expMat2)
        } else {
          sig_geneMod[[k]][[cellTypes[l]]] = list(cellType=cellTypes[l],rds_exons_t=expMat2)
        }
        
      }
    }
    
    modified_sig_geneMod[[f]] = sig_geneMod
    
  }
  
  return(modified_sig_geneMod)
  
} # End mod_sig_gM() function

#-----------------------------------------------------------#
# Step 3
# Pure Cell Type Parameter Estimation                       
#-----------------------------------------------------------#

pure_estimation = function(modified_sig_geneMod, cellTypes){
  
  pure_est = list()
  
  for(j in 1:length(modified_sig_geneMod)){
    
    if(j == 1){
      
      sig_geneMod = modified_sig_geneMod[[j]]
      
      # Only need to calcuate the pure sample parameters once
      
      sim.out = sig_geneMod[which(names(sig_geneMod)!="Sample_Info")]
      
      # Clusters with single isoforms: 
      # EXCLUDE THEM FOR THE MOMENT!
      dim_mat = matrix(0,nrow=length(sim.out),ncol=2)
      excl_clust = c()
      excl_clust2 = c()
      
      for(i in 1:(length(sim.out))){ 
        dim_mat[i,] = dim(sim.out[[i]][["X"]])
        
        if(all(dim_mat[i,]==c(1,1))){
          excl_clust = c(excl_clust,i)
        }
        if(dim_mat[i,2] == 1){
          excl_clust2 = c(excl_clust2,i)
        }
      }
      
      excl_clust_union = union(excl_clust,excl_clust2)
      if(length(excl_clust_union)>0){
        sim.new = sim.out[-c(excl_clust_union)]
      } else {
        sim.new = sim.out
      }
      
      # Optimize the Pure Sample Functions:
      ## See R/Production_Functions_PureSamp.R for Pure.apply.fun() code
      tmp.data = Pure.apply.fun(data.list = sim.new, cellTypes = cellTypes, corr_co = 1)
      
      pure_est[[j]] = tmp.data
      
    }else{
      
      tmp_j = modified_sig_geneMod[[j]][which(names(modified_sig_geneMod[[j]])!="Sample_Info")]
      
      tmp.data = pure_est[[1]]
      
      # Share calculated pure sample paramters with other list elements associated with mixtures samples
      for(clust in names(tmp.data)){
        tmp_j[[clust]][["X.fin"]] = tmp.data[[clust]][["X.fin"]]
        tmp_j[[clust]][["X.prime"]] = tmp.data[[clust]][["X.prime"]]
        tmp_j[[clust]][["l_tilde"]] = tmp.data[[clust]][["l_tilde"]]
        tmp_j[[clust]][["I"]] = tmp.data[[clust]][["I"]]
        tmp_j[[clust]][["E"]] = tmp.data[[clust]][["E"]]
        tmp_j[[clust]][["alpha.est"]] = tmp.data[[clust]][["alpha.est"]]
        tmp_j[[clust]][["beta.est"]] = tmp.data[[clust]][["beta.est"]]
        for(ct in cellTypes){
          tmp_j[[clust]][[ct]][["tau.hat"]] = tmp.data[[clust]][[ct]][["tau.hat"]]
          tmp_j[[clust]][[ct]][["gamma.hat"]] = tmp.data[[clust]][[ct]][["tau.hat"]]
        }
        
      }
      
      pure_est[[j]] = tmp_j
      
    }
    
  }
  
  return(pure_est)
  
} # End pure_estimation() function

# Downloads all count text files, computes total counts for each file
#' @export
comp_total_cts = function(directory, countData){
  
  counts_list = list()
  total_cts = numeric(length(countData))
  
  for(i in 1:length(countData)){
    countsi = read.table(file = countData[i], as.is = T)
    if(ncol(countsi) != 2){
      if(ncol(countsi) == 1 & !is.null(rownames(countsi))){
        # Perform some checks 
        if(!all(is.numeric(countsi[,1]))){
          stop(countData[i], " should have 2 columns: count and exons \n", 
               "OR should have column of counts with rownames of exons \n")
        }
        countsi[,2] = rownames(countsi)
        colnames(countsi) = NULL
        rownames(countsi) = NULL
      }else{
        stop(countData[i], " should have 2 columns: count and exons \n", 
             "OR should have column of counts with rownames of exons \n")
      }
    }
    if(all(is.numeric(countsi[,1]))){
      colNames = c("count","exons")
    }else if(all(is.numeric(countsi[,2]))){
      colNames = c("exons","count")
    }
    
    colnames(countsi) = colNames
    
    counts_list[[i]] = countsi
    counts_col = countsi[,1]
    total_cts[i] = sum(counts_col)
  }
  
  return(list(total_cts = total_cts, counts_list = counts_list))
} # End comp_total_cts() function


# export
# simControl = function(sim = FALSE, sim_num = 4, seed = NULL){
#   structure(list(sim = sim, sim_num = sim_num, seed = seed), 
#             class = c("simControl"))
# }

#' @export
optimControl = function(simple.Init = FALSE, optimType = c("nlminb")){
  structure(list(simple.Init = simple.Init, optimType = optimType),
            class = c("optimControl"))
}
