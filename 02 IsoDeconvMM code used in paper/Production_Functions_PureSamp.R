#---------------------------------------------------------------------------------#
# PURE SAMPLE: Model 2 - Original Parametrization (Logarithmic Adaptive Barrier)  #
#---------------------------------------------------------------------------------#
# whichMat: find row column indices for a feature
whichMat<-function(mat,which_idx){
  # Matrix Values
  Matr = nrow(mat)
  Matc = ncol(mat)
  
  # Transform which idx into coordinates
  col_idx = ceiling(which_idx/Matr)
  row_idx = which_idx%%Matr
  r2rep = which(row_idx==0)
  row_idx[r2rep] = Matr
  
  # for our specialized application
  tmp_mat = cbind(row_idx,col_idx)
  tmp_mat = cbind(tmp_mat,0)
  tmp_mat[,3] = c(mat[which_idx])
  idx2exc = which(apply(X=tmp_mat,MARGIN = 1,FUN = function(x){return(x[1]==x[2])}))
  
  if(length(idx2exc>0)){
    tmp_mat = tmp_mat[-c(idx2exc),]
    tmp_mat = tmp_mat[order(tmp_mat[,1]),,drop=FALSE]
  }
  
  return(tmp_mat)
}

# Picking which isoforms to exclude based on correlatoin cutoff:
Iso2Excl<-function(pairs){
  delete_iso = c()
  j=0
  while(dim(pairs)[1]>0){
    delete_iso = c(delete_iso,pairs[1,2])
    j = j+1
    
    if(dim(pairs)[1]>1){
      pairs = pairs[-c(1),,drop=FALSE]
      id2rm = which(apply(X = pairs,MARGIN = 1,FUN = function(x,del_num){return(any(x==del_num))},del_num=delete_iso[j]))
      if(length(id2rm>0)){
        pairs = pairs[-c(id2rm),,drop=FALSE]
      }
    } else {
      break
    }
  }
  return(delete_iso)
}

#----------------- Old Validated Functions--------------------------#
#--- Beta Functions ---#
ldd_beta_inv<-function(beta.curr,dim){
  JJt = matrix(1,nrow=2,ncol=2)
  A_inv = -diag(c(1/trigamma(beta.curr)))
  
  out = (1/dim)*(A_inv-(A_inv%*%JJt%*%A_inv)/((1/trigamma(sum(beta.curr)))-sum(1/trigamma(beta.curr))))
  return(out)
}

ldd_beta <- function(beta.curr,dim){
  JJt = matrix(1,nrow=2,ncol=2)
  out = dim*(trigamma(sum(beta.curr))*JJt-diag(drop(trigamma(beta.curr))))
  return(out)
}

ld_beta<-function(beta.curr,fix.tau){
  J = matrix(1,nrow=2)
  out = ncol(fix.tau)*(digamma(sum(beta.curr))*J-digamma(beta.curr))+rowSums(log(fix.tau))
  return(out)
}

#--- Alpha Functions ---#
ldd_alpha_inv<-function(alpha.curr,dim){
  JJt = matrix(1,nrow=length(alpha.curr),ncol=length(alpha.curr))
  A_inv = -diag(c(1/trigamma(alpha.curr)))
  
  out = (1/dim)*(A_inv-(A_inv%*%JJt%*%A_inv)/((1/trigamma(sum(alpha.curr)))-sum(1/trigamma(alpha.curr))))
  return(out)
}

ldd_alpha<-function(alpha.curr,dim){
  JJt = matrix(1,nrow=length(alpha.curr),ncol=length(alpha.curr))
  out = (dim)*(trigamma(sum(alpha.curr))*JJt-diag(drop(trigamma(alpha.curr))))
  return(out)
}

ld_alpha<-function(alpha.curr,fix.iso){
  J = matrix(1,nrow=nrow(fix.iso))
  out = ncol(fix.iso)*(digamma(sum(alpha.curr))*J-digamma(alpha.curr))+rowSums(log(fix.iso))
  return(out)
}

#----------------- UPDATE THE BETA VALUES --------------------------#
update.bi.v3<-function(tau.fix,init.beta){
  #--- Extract Data    ---#
  nk = ncol(tau.fix)
  
  # Update Objective/Derivative Functions:
  o.b.cfun<-function(z){
    if(any(z<=0)){
      out.lik = -Inf
    } else {
      out.lik = nk*(lgamma(sum(z))-sum(lgamma(z)))+sum((z-1)*rowSums(log(tau.fix)))
    }
    return(-out.lik)
  }
  
  g.b.cfun<-function(z){
    return(-ld_beta(beta.curr = z,fix.tau = tau.fix))
  }
  
  h.b.cfun<-function(z){
    return(-ldd_beta(beta.curr = z,dim = nk))
  }
  
  # Update Beta Values:
  auglagOut = tryCatch(nlminb(start = init.beta,objective = o.b.cfun,
                              gradient = g.b.cfun,hessian = h.b.cfun,
                              lower = 0),
                       error = function(e){
                                return(list(par=rep(NA,length(init.beta)),
                                            WARN=TRUE,convergence=100))
                              })
  
  if(auglagOut$convergence==0){
    new.beta = auglagOut$par
  } else {
    auglagOut2 = tryCatch(nlminb(start = init.beta,objective = o.b.cfun,
                               lower = 0),
                        error = function(e){
                          return(list(par=rep(NA,length(init.beta)),
                                      WARN=TRUE,convergence=100))
                        })
    
    if(auglagOut2$convergence==100){
      new.beta = rep(NA,length(init.beta))
    } else {
      new.beta = auglagOut2$par
    }
  }
  
  #------ STORE OUTPUT ------#
  return(new.beta)
}


#--------------- UPDATE THE ALPHA VALUES ----------------------#
# Revised Fit Procedure
update.ai.v3<-function(fix.iso,alpha.init){
  #------ Extracting Data Values ------#
  # Gamma Related:
  nk = ncol(fix.iso)
  
  #---------------- Update the Alpha ------------------------#
  #--- Update Objective/Derivative Functions ---#
  o.a.cfun<-function(z){
    if(any(z<=0)){
      out.lik = -Inf
    } else {
      out.lik = nk*(lgamma(sum(z))-sum(lgamma(z)))+sum((z-1)*rowSums(log(fix.iso)))
    }
    return(-out.lik)
  }
  
  g.a.cfun<-function(z){
    return(-ld_alpha(alpha.curr = z,fix.iso = fix.iso))
  }
  
  h.a.cfun<-function(z){
    return(-ldd_alpha(alpha.curr=z,dim=nk))
  }
  
  #--- Update Alpha Values ---#
  auglagOut = tryCatch(nlminb(start = alpha.init,objective = o.a.cfun,
                              gradient = g.a.cfun,hessian = h.a.cfun,lower=0),
                       error=function(e){
                              return(list(par=rep(NA,length(alpha.init)),
                                          WARN=TRUE,convergence=100))
                             })
  
  if(auglagOut$convergence==0){
    new.alpha = auglagOut$par
  } else {
    auglagOut2 = tryCatch(nlminb(start = alpha.init,
                                 objective = o.a.cfun,
                                 lower=0),
                          error=function(e){
                                  return(list(par=rep(NA,length(alpha.init)),
                                              WARN=TRUE,convergence=100))
                          })    
    if(auglagOut2$convergence==100){
      new.alpha = rep(NA,length(alpha.init))
    } else {
      new.alpha = auglagOut2$par
    }
  }
  
  #--------------- STORE THE OUTPUT ---------------------------#
  return(new.alpha)
}


#--------------- UPDATE THE GAMMA VALUES ---------------------------#
update.gamma.v.abo<-function(gamma.init,X.fin,X.prime,l_tilde,I,
                             rds.k,Xorig){
  #----------------------------------------------------------#
  # Initializing The optimization
  #----------------------------------------------------------#
  # The Method below has proven better in many cases, but much
  # more susceptible to bad points. Perhaps by improving the 
  # starting point, we improve the estimate. To do this, 
  # I introduce the constrOptim version first, then perform the
  # next step. 
  X.prime.0 = t(t(X.prime)/l_tilde[-I])
  X.fin.0   = cbind(X.prime.0,X.fin[,I])
  
  o.g.cfun<-function(z){
    gamma.probs = c(z,1-sum(z))
    if(any(gamma.probs>1)||any(gamma.probs<0)){
      out.lik = -Inf
    } else {
      out.lik = sum(rds.k[-1]*log(X.fin.0%*%matrix(c(z,1),ncol=1)))
    }
    return(out.lik)
  }
  
  g.g.cfun<-function(z){
    mu.k = X.fin.0%*%matrix(c(z,1),ncol=1)
    return(t(X.prime.0)%*%matrix(c(rds.k[-1]/mu.k),ncol=1))
  }
  
  # Update gamma Values:
  ui = rbind(diag(1,(I-1)),rep(-1,I-1))
  ci = c(rep(0,(I-1)),-1)
  new.gamma.out = tryCatch(constrOptim(theta = gamma.init,f = o.g.cfun,grad = g.g.cfun,
                                       ui = ui,ci = ci,control = list(fnscale=-1,maxit=500),
                                       method="BFGS",hessian=TRUE),error=function(e){
                                         return(list(par=c(rep((1/I),(I-1))/l_tilde[-I]),convergence=-10))
                                       })
  
  #----------------------------------------------------------#
  #
  #----------------------------------------------------------#
  o.rpg.cfun<-function(z){
    gamma.probs = c(l_tilde[-I]*exp(-z),1-sum(l_tilde[-I]*exp(-z)))
    if(any(gamma.probs>1)||any(gamma.probs<0)){
      out.lik = -Inf
    } else {
      out.lik = sum(rds.k[-1]*log(X.fin%*%matrix(c(exp(-z),1),ncol=1)))
    }
    return(out.lik)
  }
  
  g.rpg.cfun<-function(z){
    mu.k = X.fin%*%matrix(c(exp(-z),1),ncol=1)
    return((-exp(-z)*t(X.prime))%*%matrix(c(rds.k[-1]/mu.k),ncol=1))
  }
  
  hin.cfun<-function(z){
    return(-sum(l_tilde[-I]*exp(-z))+1-0.0000005)
  }
  
  hin.jac.cfun<-function(z){
    return(t(l_tilde[-I]*exp(-z)))
  }
  
  init.rp.gamma = -log(new.gamma.out$par/l_tilde[-I])
  # Not right, need to use an nlOpt function here with some new optimization routines
  # that I don't know much about. (L-BFGS is described as being similar to Newton's Method)
  new.gamma.rpg = tryCatch(auglag(par = init.rp.gamma,fn = o.rpg.cfun,gr = g.rpg.cfun,hin = hin.cfun,hin.jac = hin.jac.cfun,
                         control.optim=list(fnscale=-1),control.outer=list(trace=FALSE)),error=function(e){
                           return(list(par=-log(rep((1/I),(I-1))/l_tilde[-I]),convergence=-10))
                         })
  
  new.gamma = exp(-new.gamma.rpg$par)
  new.gamma.err = new.gamma.rpg$convergence
  
  return(list(new.gamma=new.gamma,gamma.err = new.gamma.err))
}

#--------------- OPTIMIZE IN THE PURE SAMPLE ------------------------#
update.pure.alpha.beta.sing<-function(p.data,X.fin,X.prime,l_tilde,I,Xorig){
  #----- (0): Extract Data ----#
  rds.k = p.data[["rds_exons"]]
  
  #----- (1): Tau Updates -----#
  if(ncol(rds.k) > 1){
    tau.hat = colSums(rds.k[-1,])/colSums(rds.k)
  }else{
    # warn_msg = sprintf("Number of reference pure samples for cell type %s is 1 \n IsoDeconvMM suggests using 2 or more reference samples per cell type \n", p.data$cellType)
    # warning(warn_msg, immediate. = T)
    tau.hat = sum(rds.k[-1,])/sum(rds.k)
  }
  
  
  #--- (2): Initialize Beta ---#
  mu.1 = mean(tau.hat)
  mu.2 = mean(tau.hat^2)
  b1.init = (mu.1^2)*(1-mu.1)/(mu.2-mu.1^2)-mu.1
  b2.init = ((1-mu.1)/mu.1)*b1.init
  b.j.init = matrix(c(b1.init,b2.init),nrow=2)
  fix.tau = matrix(c(tau.hat,1-tau.hat),nrow=2,byrow=TRUE)
  
  #--- (3): Update the Beta ---#
  beta.hat = update.bi.v3(tau.fix=fix.tau,
                          init.beta = b.j.init)
  
  #--- (4): Gamma Updates ---#
  gamma.init = rep((1/I),(I-1))/l_tilde[-I]
  gamma.hat.t = apply(X = rds.k,MARGIN = 2,FUN = update.gamma.v.abo,
                      X.fin=X.fin,X.prime=X.prime,l_tilde=l_tilde,I=I,
                      gamma.init=gamma.init,Xorig=Xorig)
  gamma.hat.vals = matrix(0,ncol=ncol(rds.k),nrow=(I-1))
  gamma.hat.conv = matrix(0,ncol=ncol(rds.k),nrow=1)
  for(j in 1:ncol(rds.k)){
    gamma.hat.vals[,j] = gamma.hat.t[[j]]$new.gamma
    gamma.hat.conv[,j] = gamma.hat.t[[j]]$gamma.err
  }
  fix.iso = rbind(l_tilde[-I]*gamma.hat.vals,1-colSums(l_tilde[-I]*gamma.hat.vals))
  
  # Excluding Bad Samples:
  cols2rm = which(gamma.hat.conv==-10)
  if(length(cols2rm)==0){
    fix.iso.use = fix.iso
  } else {
    fix.iso.use = fix.iso[,-c(cols2rm)]
  }
  
  # Allow for slight error in clusters where it is estimated
  # that only a single isoform is used
  fix.iso.use = apply(X=fix.iso.use,MARGIN=2,FUN = function(z){
    if(any(z>=0.98)){
      idx = which(z>=0.98)
      z[idx] = 0.98
      z[-idx] = 0.02/(length(z)-1)
    } else {
      
    }
    return(z)
  })
  
  #--- (5): Initialize Alpha ---#
  mu = apply(X = fix.iso.use,MARGIN = 1,FUN = mean)
  mu.2 = apply(X = (fix.iso.use^2),MARGIN = 1,FUN = mean)
  if(all((mu.2-mu^2)<1e-8)){
    # Little Variation is observed in our data
    # adjust alpha penalty to allow for this
    # Value chosen so that at most there is 1% SD on 
    # isoform parameters
    alpha.hat = mu*(2500)
  } else {
    a.sum.init = ((mu*(1-mu))/(mu.2-mu^2))-1
    a.j.init   = a.sum.init*mu
    a.j.init[which(a.j.init==Inf)] = 1e5
    
    #--- (6): Update the Alpha ---#
    alpha.hat = update.ai.v3(fix.iso = fix.iso.use,
                             alpha.init = a.j.init)
  }
  
  return(list(tau.hat=tau.hat,beta.hat=beta.hat,gamma.hat=gamma.hat.vals,alpha.hat=alpha.hat))
}

update.pure.alpha.beta.all<-function(data.list,cellTypes,X.fin,X.prime,l_tilde,I){
  data.list.pure = data.list[cellTypes]
  
  out.pure.update = lapply(X = data.list.pure,FUN = update.pure.alpha.beta.sing,
                           X.fin=X.fin,X.prime=X.prime,l_tilde=l_tilde,I=I,Xorig=data.list$X)
  return(out.pure.update)
}

#----------------- Function for Applying to each Cluster --------------------------------#
Pure.update.cycle<-function(z,cellTypes,corr_co){
  #----------------------------------------------------------#
  # STEP 0: Prepare X-matrices                               #
  #----------------------------------------------------------#
  X = z$X
  I = ncol(X)
  E = nrow(X)
  l_tilde = colSums(X)
  
  #------------ Removing Correlated Isoforms ----------------#
  IsoCorr = cor(X)
  
  iso2d = NULL
  
  if(any(IsoCorr>corr_co)){
    IsoCorr[lower.tri(x=IsoCorr,diag=TRUE)]=0
    pairs_v = whichMat(mat = IsoCorr,which_idx = which(IsoCorr>corr_co))
    iso2d = Iso2Excl(pairs = pairs_v)
  }
  
  if(length(iso2d)>0){
    X = X[,-c(iso2d)]
  }
  
  l_tilde = colSums(X)
  I = ncol(X)
  
  #-------------------------------------------------------------------------#
  X.fin = X[,-I] - dupcol(mult_vec = l_tilde[-I]/l_tilde[I],col_vec = X[,I])
  X.prime = X.fin
  X.fin = cbind(X.fin,X[,I]/l_tilde[I])
  
  z[["X.fin"]] = X.fin
  z[["X.prime"]] = X.prime
  z[["l_tilde"]] = l_tilde
  z[["I"]] = I
  z[["E"]] = E
  
  #----------------------------------------------------------#
  # STEP 1: Set Up List Necessities                          #
  #----------------------------------------------------------#
  z[["alpha.est"]] = matrix(0,ncol=length(cellTypes),nrow=I)
  colnames(z[["alpha.est"]]) = cellTypes
  
  z[["beta.est"]] = matrix(0,ncol=length(cellTypes),nrow=2)
  colnames(z[["beta.est"]]) = cellTypes
  
  #---------------- Update Pure Sample Values -------------------#
  out.final = update.pure.alpha.beta.all(data.list = z,cellTypes = cellTypes,
                                         X.fin = X.fin,X.prime = X.prime,
                                         l_tilde = l_tilde,I = I)
  
  # Update the list:
  for(i in 1:length(cellTypes)){
    z[[cellTypes[i]]]$tau.hat = out.final[[cellTypes[i]]]$tau.hat
    z[[cellTypes[i]]]$gamma.hat = out.final[[cellTypes[i]]]$gamma.hat
    z[["alpha.est"]][,cellTypes[i]] = out.final[[cellTypes[i]]]$alpha.hat
    z[["beta.est"]][,cellTypes[i]] = out.final[[cellTypes[i]]]$beta.hat
  }
  
  return(z)
}

#' @import alabama
#' @import gtools
Pure.apply.fun<-function(data.list,cellTypes,corr_co){
  out.fin = lapply(X = data.list,FUN = Pure.update.cycle,cellTypes=cellTypes,corr_co=corr_co)
  return(out.fin)
}

dupcol<-function(mult_vec,col_vec){
  #---------------- VARIABLE DESCRIPTION ----------------#
  # *col_vec  : column vector to construct matrix
  # *Mult_Vec : Vector of values used to multiply col_vec
  l1 = length(mult_vec)
  l2 = length(col_vec)
  out.mat = t(t(matrix(col_vec,nrow=l2,ncol=l1))*mult_vec)
  return(out.mat)
}