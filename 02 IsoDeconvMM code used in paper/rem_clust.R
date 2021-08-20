
rem_clust<-function(geneMod,co,min_ind){
  indices = which(names(geneMod)!="Sample_Info")
  geneMod_labels = paste("y_",geneMod[["Sample_Info"]]$Label,sep="")
  co_ct = c()
  
  for(i in indices){
    tot_all = rep(0,length(geneMod_labels))
    for(j in 1:length(geneMod_labels)){
      tot_all[j] = sum(geneMod[[i]][[geneMod_labels[j]]])
    }
    if(min_ind==1 && min(tot_all)<co){co_ct=c(co_ct,i)}
    if(min_ind==0 && median(tot_all)<co){co_ct=c(co_ct,i)}
    if(all(dim(as.matrix(geneMod[[i]][["X"]]))==c(1,1))==1){co_ct=c(co_ct,i)}
    if(any(dim(as.matrix(geneMod[[i]][["X"]]))==0)){co_ct = c(co_ct,i)}
  }
  co_ct = unique(co_ct)
  
  if(length(co_ct)!=0){geneMod = geneMod[-co_ct]}
  
  return(geneMod)
}