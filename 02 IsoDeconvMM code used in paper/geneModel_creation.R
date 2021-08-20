
geneModel_creation = function(concat_geneMod, fragSizeFile, cellTypes, labels,
                            knownIsoforms, readLen, lmax, eLenMin, total_cts){
  
  # Define variables
  isoAll = knownIsoforms
  sam_names = labels
  isoAll = knownIsoforms
  
  #-------------------------------------------------------------------------------------------------------------------------------------------#
  # GENERATING PDDIST:                                                                                                                        #
  #-------------------------------------------------------------------------------------------------------------------------------------------#
  # Generates an estimate of the distribution of fragment lengths necessary for the computation of effective length.
  
  pdDist_gen <- function(fragSizeFile,lmax){
    
    md = fragSizeFile
    
    pd = rep(0, lmax)
    w2 = which(md$Len <= lmax)
    pd[md$Len[w2]] = md$Freq[w2]
    pdDist = pd/sum(pd)
    return(pdDist)
  }
  
  pdDist = pdDist_gen(fragSizeFile,lmax)
  
  #-------------------------------------------------------------------------------------------------------------------------------------------#
  # CALL TO GENEMODEL:                                                                                                                        #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Calls to the geneModel function for each transcript cluster.                                                                              #
  #-------------------------------------------------------------------------------------------------------------------------------------------#
  
  fin_geneMod = list()
  nms = names(concat_geneMod)
  sep = "------"
  mix_sams = sam_names[which(tolower(cellTypes)=="mix")]
  
  for (i in 1:length(nms)) {
    if (i%%100 == 0) {
      message(sep, i, "  ", date(), sep)
    }
    nm1 = nms[i]
    ge1 = concat_geneMod[[nm1]]
    isoforms = isoAll[[nm1]]
    
    # Call to geneModel function: ge1 is a list with components $count and $info for a single transcript
    # cluster. isoforms is a matrix which details all of the isoforms used in a transcript cluster with
    # indicators for whether a particular exon is used by an isoform.
    
    # See "R/geneModel_multcell_edit.R" file for geneModel_multcell_Edit() function code
    gm1 = geneModel_multcell_Edit(ge1, d = readLen, pdDist, isoforms, lmax, 
                                  eLenMin, verbose=1,sam_names=sam_names,mix_sams=mix_sams)
    fin_geneMod[[nm1]] = gm1
  }
  
  
  info_mat = data.frame(Label = labels, Cell_Type = cellTypes, Total = total_cts, stringsAsFactors = FALSE)
  
  fin_geneMod[["Sample_Info"]] = info_mat
  
  return(fin_geneMod)
  
}