#-----------------------------------------------------------------------#
# PRODUCTION FUNCTIONS                                                  #
#-----------------------------------------------------------------------#
# NAME:                                                                 #
#   Production Functions.R                                              #
# DATE:                                                                 #
#   12/30/2015                                                          #
# VERSION:                                                              #
#   R.2.13.5                                                            #
# PROGRAMMER:                                                           #
#   Douglas Roy Wilson, Jr.                                             #
#-----------------------------------------------------------------------#
# DESCRIPTION:                                                          #
#   Compiles all of the necessary functions for the fit of the staged   #
#   model estimation approach. This document compiles and organizes the #
#   fit in the mixture tissue only.                                     #
#-----------------------------------------------------------------------#

#-----------------------------------------------------------------------#
# SUPPORTING FUNCTIONS                                                  #
#-----------------------------------------------------------------------#

#------------------------------------------------------------------#
# Mixture Likelihood for Propn. and Taus                           #
#------------------------------------------------------------------#
mix.full.stage.lik.pt<-function(beta,pm.rds,p,tau){
  #====================================================================================#
  # - beta    : matrix of beta hyperparameters for each cell type                      #
  #------------------------------------------------------------------------------------#
  # - pm.rds  : vector of posterior mean total cell type counts (of length nk)         #
  #------------------------------------------------------------------------------------#
  # - p       : vector of length nk specifying proportions of each cell type           #
  #------------------------------------------------------------------------------------#
  # - tau     : vector of length nk specifying tau parameters for each cell type       #
  #------------------------------------------------------------------------------------#
  # - gamma   : Matrix (I-1 x nk) of the sample specific isoform expression parameters #
  #====================================================================================#
  
  #--- Negativity Criteria Check ---#
  # relaxes problems at boundary to avoid optimization issues
  # 0.000001 represents a negligible proportion in your sample
  # and could be set to 0 without any harm being done.
  cond.p = any(p<(-1e-6))||any(p>(1+1e-6))
  cond.tau = any(tau<0)||any(tau>1)
  
  if(any(cond.p,cond.tau)){
    out.lik = -Inf
    return(out.lik)
  } else if(any(p<0)||any(p>1)){
    idx0 = which(p<0)
    idx1 = which(p>1)
    p[idx0] = 1e-16
    p[idx1] = 1-(1e-16)
  } 
  
  #--- Multinomial Component ---#
  trunc.prob = c(p*tau/sum(p*tau))
  m.lik.o = pm.rds*log(trunc.prob)
  
  #--- Tau Dirichlet Component ---#
  tau.lik.o = lgamma(colSums(beta))-colSums(lgamma(beta))+colSums((beta-1)*log(matrix(c(tau,1-tau),nrow=2,byrow = TRUE)))
  
  #--- Complete Value ---#
  out.lik = sum(m.lik.o)+sum(tau.lik.o)
  
  #--- Return Output ---#
  return(out.lik)
}

ActLik4test<-function(p,tau,gamma.est,cts,X.fin,alpha.est,beta.est,l_tilde){
  # Slight Fudging to allow for comparison:
  iso.probs = gamma.est*l_tilde
  iso.probs = rbind(iso.probs,1-colSums(iso.probs))
  
  for(i in 1:ncol(alpha.est)){
    iso.ind = any(iso.probs[,i]==0)
    if(iso.ind==1){
      w2change = which(iso.probs[,i]==0)
      iso.probs[w2change,i] = 0.0025
    }
  }
  
  if(any(p< (-1e-6))||any(p>(1+1e-6))){
    out.lik=-Inf
    return(out.lik)
  } 
  
  if(any(p<0)||any(p>1)){
    idx0 = which(p<0)
    idx1 = which(p>1)
    p[idx0] = 0.0025
    p[idx1] = 1-0.0025
  } 
  
  gamma.est.f = iso.probs[-c(nrow(iso.probs)),]
  gamma.est.f = gamma.est.f/l_tilde
  
  condProbs = X.fin%*%rbind(gamma.est,1)
  condProbs = t(t(condProbs)*c(p*tau))/sum(p*tau)
  condProbsfin = rowSums(condProbs)
  
  out.m.lik = cts[-1]*log(condProbsfin)
  out.a.lik = lgamma(colSums(alpha.est))-colSums(lgamma(alpha.est))+colSums((alpha.est-1)*log(iso.probs))
  out.b.lik = lgamma(colSums(beta.est))-colSums(lgamma(beta.est))+colSums((beta.est-1)*log(matrix(c(tau,1-tau),nrow=2,byrow=TRUE)))
  
  out.lik = sum(out.m.lik)+sum(out.a.lik)+sum(out.b.lik)
  return(out.lik)
}

applyLik<-function(clusterData){
  p = clusterData$mix$p.est
  tau = clusterData$mix$tau.est
  gamma.est = matrix(clusterData$mix$gamma.est,nrow=(clusterData$I-1),ncol=length(tau))
  cts = clusterData$mix$rds_exons_t
  X.fin = clusterData$X.fin
  alpha.est = clusterData$alpha.est
  beta.est  = clusterData$beta.est
  l_tilde   = clusterData$l_tilde[-c(clusterData$I)]
  
  return(ActLik4test(p = p,tau = tau,gamma.est = gamma.est,cts = cts,X.fin = X.fin,alpha.est = alpha.est,beta.est = beta.est,l_tilde = l_tilde))
}

#------------------------------------------------------------------#
# Compute Posterior Means                                          #
#------------------------------------------------------------------#
compute.pmeans<-function(X_fin,mix_list){
  # Inside Gene of Interest:
  condprobs = X_fin%*%rbind(mix_list$gamma.est,1)
  condprobs = t(t(condprobs)*c(mix_list$p.est*mix_list$tau.est/sum(mix_list$p.est*mix_list$tau.est)))
  pm.I = (condprobs/rowSums(condprobs))*c(mix_list$rds_exons_t[-1])
  
  # Output:
  return(pm.I)
}

#-----------------------------------------------------------------------#
# FITTING Functions                                                     #
#-----------------------------------------------------------------------#
#------------------------------------------------------------------#
# Update the Fit of the proportions and tau values                 #
#------------------------------------------------------------------#
# Originally considered for the log adaptive barrier fitting, the data
# itself contains natural barriers to the taus and proportions falling
# outside of the desired ranges. Thus, for preliminary evaluation, I 
# just allow R code to perform a Quasi-Newton optimization method that
# approximates Newton's Method. Approximates Hessian without computing.

lab_fit_update_v2<-function(cdata,nk){
  #--- Extract Cell Type Specific Quantities ---#
  X.fin = cdata$X.fin
  X.prime = cdata$X.prime
  I = cdata$I
  l_tilde = cdata$l_tilde
  
  pm.rds = colSums(cdata[["mix"]][["pm.rds.exons"]])
  
  alpha.mat = cdata[["alpha.est"]]
  beta.mat = cdata[["beta.est"]]
  
  #--- Define the Likelihood function ---#
  lik_clust_func<-function(z){
    p.tmp = z[c(1:(nk-1))]
    p.tmp = c(p.tmp,1-sum(p.tmp))
    tau.tmp = exp(-z[-c(1:(nk-1))])
    
    out.param = mix.full.stage.lik.pt(beta = beta.mat,
                                      pm.rds = pm.rds,
                                      p = p.tmp,tau = tau.tmp)
    return(-out.param)
  }
  
  #--- Define the Cluster Gradient ---#
  grad_clust_func<-function(z){
    p.tmp = z[c(1:(nk-1))]
    p.f.tmp = c(p.tmp,1-sum(p.tmp))
    tau.tmp = exp(-z[-c(1:(nk-1))])
    
    #--- Prop. Grad. ---#
    p.grad = (pm.rds[-nk]/p.tmp)-rep(sum(pm.rds),(nk-1))*(tau.tmp[-nk]-tau.tmp[nk])/sum(tau.tmp*p.f.tmp)-
      rep((pm.rds[nk]/p.f.tmp[nk]),(nk-1))
    
    #---   Tau Grad. ---#
    t.grad = -(pm.rds+beta.mat[1,]-1)+(beta.mat[2,]-1)*(tau.tmp/(1-tau.tmp))+
      sum(pm.rds)*p.f.tmp*tau.tmp/(sum(p.f.tmp*tau.tmp))
    
    #--- Ouput Grad. ---#
    return(c(-p.grad,-t.grad))
  }
  
  hinFunc<-function(z){
    outVec = rep(0,nk)
    outVec[c(1:(nk-1))] = z[c(1:(nk-1))]-(0.0025)
    outVec[nk] = 1-(0.0025)-sum(z[c(1:(nk-1))])
    return(outVec)
  }
  
  hinJac<-function(z){
    outJac = matrix(0,nrow=nk,ncol=length(z))
    outJac[c(1:(nk-1)),c(1:(nk-1))] = diag(1,nrow=(nk-1))
    outJac[nk,c(1:(nk-1))] = -1
    return(outJac)
  }
  
  #--- Run the Optimization Routine ---#
  # Perhaps try augLag here.
  est.param = auglag(par = c(cdata[["mix"]]$p.est[-nk],-log(cdata[["mix"]][["tau.est"]])),
                     fn = lik_clust_func,gr = grad_clust_func,hin = hinFunc,hin.jac = hinJac,
                     control.optim = list(maxit=1000),control.outer=list(ilack.max=6,trace=FALSE))
  p.est = est.param$par[c(1:(nk-1))]
  t.est = exp(-est.param$par[-c(1:(nk-1))])
  
  return(list(p.est = p.est,t.est = t.est))
}

#------------------------------------------------------------------#
# UPDATE Isoform values                                            #
#------------------------------------------------------------------#

update.iso.STG2<-function(X.fin,X.prime,I,l_tilde,
                          alphas.prev.rds,optimType){
  #-----------------------------------------------------------------#
  # EXTRACT alphas and reads                                        #
  #-----------------------------------------------------------------#
  alphas.curr   = alphas.prev.rds[1:I]
  gamma.init = alphas.prev.rds[c((I+1):(2*I-1))]
  pm.rds.exons  = alphas.prev.rds[-c(1:(2*I-1))]
  
  #-----------------------------------------------------------------#
  # FUNCTION DEVELOPMENT                                            #
  #-----------------------------------------------------------------#
  o.g.CFUN<-function(z){
    iso.tmp = exp(-z)
    gamma.probs = c(l_tilde[-I]*iso.tmp,1-sum(l_tilde[-I]*iso.tmp))
    if(any(gamma.probs>(1+1e-6))||any(gamma.probs<(-1e-6))){
      out.lik = -Inf
    } else if(any(gamma.probs>1)||any(gamma.probs<0)){
      idx1 = which(gamma.probs>1)
      idx0 = which(gamma.probs<0)
      gamma.probs[idx0] = 1e-6
      gamma.probs[idx1] = 1-1e-6
      gamma.probs = gamma.probs/sum(gamma.probs)
      
      iso.tmp = gamma.probs[-I]/l_tilde[-I]
      
      out.lik = sum(pm.rds.exons*log(X.fin%*%matrix(c(iso.tmp,1),ncol=1)))+lgamma(sum(alphas.curr))-
        sum(lgamma(alphas.curr))+sum((alphas.curr-1)*log(gamma.probs))
    } else {
      out.lik = sum(pm.rds.exons*log(X.fin%*%matrix(c(iso.tmp,1),ncol=1)))+lgamma(sum(alphas.curr))-
        sum(lgamma(alphas.curr))+sum((alphas.curr-1)*log(gamma.probs))
    }
    return(-out.lik)
  }
  
  g.g.CFUN<-function(z){
    iso.tmp = exp(-z)
    
    means = X.fin%*%matrix(c(iso.tmp,1),ncol=1)
    grad.out = diag(-c(iso.tmp),(I-1))%*%t(X.prime)%*%matrix(c(pm.rds.exons/means),ncol=1)-
      matrix(c((alphas.curr[-I]-1)),ncol=1)+
      matrix(c(((alphas.curr[I]-1)/(1-sum(l_tilde[-I]*iso.tmp)))*l_tilde[-I]*iso.tmp),ncol=1)
    return(-grad.out)
  }
  
  h.g.CFUN<-function(z){
    iso.tmp = exp(-z)
    gamma.probs = l_tilde[-I]*iso.tmp
    
    tmp.mat = t(t(X.prime)*c(iso.tmp))
    
    means = X.fin%*%matrix(c(iso.tmp,1),ncol=1)
    hess.out = -((alphas.curr[I]-1)/((1-sum(gamma.probs))^2))*matrix(gamma.probs,ncol=1)%*%matrix(gamma.probs,nrow=1)-
      ((alphas.curr[I]-1)/(1-sum(gamma.probs)))*diag(c(gamma.probs),(I-1))-t(tmp.mat)%*%diag(c(pm.rds.exons/(means^2)))%*%(tmp.mat)+
      diag(c(colSums(tmp.mat*c(pm.rds.exons/means))),nrow=(I-1))
    return(-hess.out)
  }
  
  if(optimType=="nlminb"){
    out.iso.r = nlminb(start = -log(gamma.init),objective = o.g.CFUN,gradient = g.g.CFUN,hessian = h.g.CFUN)
    # out.iso.r = nlminb(start = -log(gamma.init),objective = o.g.CFUN,gradient = g.g.CFUN,hessian = h.g.CFUN, lower = ?)
  } else if(optimType=="BFGS"){
    out.iso.r = optim(par = -log(gamma.init),fn = o.g.CFUN,gr = g.g.CFUN,method = "BFGS",control = list(maxit=5000))
  } else if(optimType=="NM"){
    out.iso.r = optim(par = -log(gamma.init),fn = o.g.CFUN,gr = g.g.CFUN,control = list(maxit=5000))
  } else if(optimType=="BFGSWB"){
    hin.cfun<-function(z){
      return(-sum(l_tilde[-I]*exp(-z))+0.99975)
    }
    
    hin.jac.cfun<-function(z){
      return(t(l_tilde[-I]*exp(-z)))
    }
    
    # out.iso.r = tryCatch(auglag(par = -log(gamma.init),fn = o.g.CFUN,gr = g.g.CFUN,hin = hin.cfun,hin.jac = hin.jac.cfun,
    #                                 control.optim=list(fnscale=1),control.outer=list(trace=FALSE)),error=function(e){
    #                                   return(list(par=-log(rep((1/I),(I-1))/l_tilde[-I]),convergence=-10))
    #                                 })
    out.iso.r = auglag(par = -log(gamma.init),fn = o.g.CFUN,gr = g.g.CFUN,hin = hin.cfun,hin.jac = hin.jac.cfun,
                                control.optim=list(fnscale=1),control.outer=list(trace=FALSE))
  }
  
  # Isoform Revision
  isoprobs = l_tilde[-I]*exp(-out.iso.r$par)
  isoprobst = c(isoprobs,1-sum(isoprobs))
  
  if(any(isoprobst<0)||any(isoprobst>1)){
    idx1 = which(isoprobst>1)
    idx0 = which(isoprobst<0)
    
    isoprobst[idx0] = 0.00025
    isoprobst[idx1] = 0.99975
    
    isoprobst = isoprobst/sum(isoprobst)
    
    iso.out = isoprobst[-I]/l_tilde[-I]
  } else {
    iso.out = exp(-out.iso.r$par)
  }
   
  return(iso.out)
}

#-----------------------------------------------------------------------#
# WRAPPER functions                                                     #
#-----------------------------------------------------------------------#
# Let's just set up for the two cell type case. If length celltypes exceeds 2,
# let's output a warning.  

applyMI<-function(x,cdata,cellTypes,p.co,t.co,g.co,optimType,iter.cut){
  tmp.Data = STG.Update_Cluster.Single(cdata = cdata,cellTypes = cellTypes,
                                       p.co = p.co,t.co =t.co ,g.co = g.co,
                                       optimType = optimType,test.init = c(x),
                                       iter.cut = iter.cut)
  if(tmp.Data$WARN%in%c(4,5)){
    out.Lik = -Inf
  } else {
    out.Lik = applyLik(tmp.Data)
  }
  
  return(list(out.Lik = out.Lik,p.est = tmp.Data$mix$p.est))
}

STG.Update_Cluster.SingMI<-function(cdata,cellTypes,p.co,t.co,g.co,optimType,iter.co,simple.Init,initPts){
  
  # initPts values set up in isoDeconvMM() function
  test.init = initPts
  
  
  if(simple.Init==TRUE){
    message("Applying First Round!")
    
    outLikVals_List = apply(X = test.init,MARGIN = 1,FUN = applyMI,
                       cdata=cdata,cellTypes=cellTypes,p.co=1e-4,t.co=1e-7,g.co=1e-4,
                       optimType="BFGSWB",iter.cut=7500)
    
    outLikVals = rep(0,length(outLikVals_List))
    est.Props  = matrix(0,nrow=length(outLikVals_List),ncol=length(cellTypes))
    for(i in 1:length(outLikVals_List)){
      outLikVals[i] = outLikVals_List[[i]]$out.Lik
      est.Props[i,] = outLikVals_List[[i]]$p.est
    }
    
    spt2use = which.max(outLikVals)
    
    maxDiff = (max(outLikVals,na.rm = TRUE)-min(outLikVals,na.rm = TRUE))
    
    if(is.na(maxDiff)){
      tmp.out = cdata
      tmp.out[["WARN"]] = 4
    } else if(abs(maxDiff/max(outLikVals,na.rm=TRUE))<0.001){
      rangeEst = apply(X = est.Props,MARGIN = 2,FUN = function(x){
        return(max(x,na.rm = TRUE)-min(x,na.rm=TRUE))
      })
      
      if(max(rangeEst)<0.001){
        # All proportions estimated similarly, pick a start point in the middle
        # and go with it.
        tmp.out = STG.Update_Cluster.Single(cdata=cdata,cellTypes=cellTypes,p.co=p.co,t.co=t.co,
                                            g.co=g.co,optimType=optimType,
                                            test.init=c(test.init[spt2use,]))
      } else {
        message("Seemingly Flat Likelihood! More Detail Needed to Pick Starting Point!")
        tmp.out = STG.Update_Cluster.Single(cdata=cdata,cellTypes=cellTypes,p.co=p.co,t.co=t.co,
                                            g.co=g.co,optimType=optimType,
                                            test.init=c(test.init[1,]))
        for(j in 2:nrow(test.init)){
          message("Examining start point ",j,"!")
          tmp.new = STG.Update_Cluster.Single(cdata=cdata,cellTypes=cellTypes,p.co=p.co,t.co=t.co,
                                              g.co=g.co,optimType=optimType,
                                              test.init=c(test.init[j,]))
          
          elements_NA = sapply(tmp.new, FUN = function(x) sum(is.na(x)))
          if(sum(elements_NA) > 0){
            # If any elements of tmp.new are NA, then skip
            message("NA values occurred in fit ", j, ", skipping to next start point")
            next
          }else if(sum(sapply(tmp.out, FUN = function(x) sum(is.na(x))))){
            # If first tmp.out has NA values but tmp.new is fine, replace tmp.out with tmp.new 
            # and skip to next start point
            tmp.out = tmp.new
            next
          }
          
          # If applyLik output is NA, skip
          if(is.na(applyLik(tmp.new))){
            message("NA values occurred in fit ", j, ", skipping to next start point")
            next
          }else if(is.na(applyLik(tmp.out))){ # in case tmp.out results in NA value
            tmp.out = tmp.new
            next
          }
          
          if(applyLik(tmp.new)>applyLik(tmp.out)){
            tmp.out = tmp.new
          }
        }
      }
    } else {
      tmp.out = STG.Update_Cluster.Single(cdata=cdata,cellTypes=cellTypes,p.co=p.co,t.co=t.co,
                                          g.co=g.co,optimType=optimType,
                                          test.init=c(test.init[spt2use,]))
    }
    
  } else {
    message("Simple Init Not Performed! Full Fit for each start point!")
    tmp.out = STG.Update_Cluster.Single(cdata=cdata,cellTypes=cellTypes,p.co=p.co,t.co=t.co,
                                        g.co=g.co,optimType=optimType,
                                        test.init=c(test.init[1,]))
    
    if(nrow(test.init)>1){
      for(j in 2:nrow(test.init)){
        message("Examining start point ",j,"!")
        tmp.new = STG.Update_Cluster.Single(cdata=cdata,cellTypes=cellTypes,p.co=p.co,t.co=t.co,
                                            g.co=g.co,optimType=optimType,
                                            test.init=c(test.init[j,]))
        
        
        elements_NA = sapply(tmp.new, FUN = function(x) sum(is.na(x)))
        if(sum(elements_NA) > 0){
          # If any elements of tmp.new are NA, then skip
          message("NA values occurred in fit ", j, ", skipping to next start point")
          next
        }else if(sum(sapply(tmp.out, FUN = function(x) sum(is.na(x))))){
          # If first tmp.out has NA values but tmp.new is fine, replace tmp.out with tmp.new 
          # and skip to next start point
          tmp.out = tmp.new
          next
        }
        
        # If applyLik output is NA, skip
        if(is.na(applyLik(tmp.new))){
          message("NA values occurred in fit ", j, ", skipping to next start point")
          next
        }else if(is.na(applyLik(tmp.out))){ # in case tmp.out results in NA value
          tmp.out = tmp.new
          next
        }
        
        if(applyLik(tmp.new)>applyLik(tmp.out)){
          tmp.out = tmp.new
        }
        
      }
    }
    
  }
  
  return(tmp.out)
}

STG.Update_Cluster.Single<-function(cdata,cellTypes,p.co,t.co,g.co,optimType,test.init,iter.cut=15000){
  #----------------------------------------------------#
  # EM Loop for Estimation of Proportions              #
  #----------------------------------------------------#
  # EM convergence decided by the parameters of interest
  # (the proportions). 
  p.conv = TRUE
  t.conv = TRUE
  g.conv = TRUE
  mix.tot = sum(cdata[["mix"]][["rds_exons_t"]][-1])
  iter.co = 0
  warn.clust=0
  
  #----------------------------------------------------#
  # Check to make sure that we are properly ordered    #
  #----------------------------------------------------#
  #--- Placeholders ---#
  cdata[["mix"]][["gamma.est"]] = matrix(0,nrow=(cdata$I-1),
                                         ncol=length(cellTypes))
  colnames(cdata[["mix"]][["gamma.est"]]) = cellTypes
  cdata[["mix"]][["tau.est"]] = matrix(0,nrow = 1,ncol=length(cellTypes))
  colnames(cdata[["mix"]][["tau.est"]]) = cellTypes
  
  #--- Fix Order ---#
  cdata[["alpha.est"]] = cdata[["alpha.est"]][,cellTypes]
  cdata[["beta.est"]]  = cdata[["beta.est"]][,cellTypes]
  
  #--- Initialized ---#
  #cdata[["mix"]][["p.est"]] = matrix(1/length(cellTypes),nrow=1,ncol=length(cellTypes))
  cdata[["mix"]][["p.est"]] = matrix(test.init,nrow=1,ncol=length(cellTypes))
  cdata[["mix"]][["tau.est"]] = cdata[["beta.est"]][1,]/colSums(cdata[["beta.est"]])
  cdata[["mix"]][["gamma.est"]] = t(t(cdata[["alpha.est"]])/colSums(cdata[["alpha.est"]]))[-c(cdata$I),]/cdata$l_tilde[-cdata$I]
  cdata[["mix"]][["pm.rds.exons"]] = matrix(0,nrow=(cdata$E),ncol=length(cellTypes))
  
  if(any(is.na(cdata[["alpha.est"]]))||any(is.na(cdata[["beta.est"]]))){
    warn.clust = 5
    TroubleInd=NULL
  } else {
    TroubleInd<-try(
      while(any(p.conv,t.conv,g.conv)){
        #--------------------------------------------------#
        # STEP 1 (E-Step): UPDATE the Posterior Means      #
        #--------------------------------------------------#
        cdata[["mix"]][["pm.rds.exons"]] = compute.pmeans(X_fin = cdata$X.fin,
                                                          mix_list = cdata[["mix"]])
        
        #--------------------------------------------------#
        # STEP 2.A (M-Step): UPDATE tau and proportions    #
        #--------------------------------------------------#
        tmp.tp.out = lab_fit_update_v2(cdata = cdata,
                                       nk = length(cellTypes))
        
        #--------------------------------------------------#
        # STEP 2.B (M-Step): UPDATE the Gamma values       #
        #--------------------------------------------------#
        alphas.prev.rds = rbind(cdata[["alpha.est"]],
                                cdata[["mix"]][["gamma.est"]],
                                cdata[["mix"]][["pm.rds.exons"]])
        
        new.gamma = apply(X = alphas.prev.rds,MARGIN = 2,FUN = update.iso.STG2,
                          X.fin = cdata$X.fin, X.prime = cdata$X.prime,
                          I = cdata$I, l_tilde = cdata$l_tilde,optimType = optimType)
        
        #--------------------------------------------------#
        # CHECK Convergence of the Values / STORE          #
        #--------------------------------------------------#
        #----- Inflation Factors --------------------------#
        # Computing % change causes problems if p and tau have zero entries
        # Best to add inflation factor to control for this problem. 
        
        # P threshold chosen since this gets into the 0.001%
        # In mixture sample. We don't much care about that.
        
        # Tau threshold chosen since at this level we would 
        # average less than 1 read per 10,000,000 total reads
        # in the cluster. Seems a sufficiently small bump. 1 
        # read is approximately zero in my book because of error.
        
        t.inflation = 1e-8
        
        #----- Convergence Computations -------------------#
        new.p = c(tmp.tp.out$p.est,1-sum(tmp.tp.out$p.est))
        p.diff = max(abs(new.p-cdata[["mix"]]$p.est))
        cdata[["mix"]][["p.est"]] = matrix(c(tmp.tp.out[["p.est"]],1-sum(tmp.tp.out[["p.est"]])),nrow=1,ncol=length(cellTypes))
        p.conv = (p.diff>p.co)
        
        tau.diff = max(abs(tmp.tp.out$t.est-cdata[["mix"]]$tau.est))
        cdata[["mix"]][["tau.est"]] = tmp.tp.out[["t.est"]]
        t.conv = (tau.diff>t.co)
        
        gamma.diff = max(abs(cdata[["mix"]][["gamma.est"]]-new.gamma)*cdata$l_tilde[-cdata$I])
        cdata[["mix"]][["gamma.est"]] = new.gamma
        g.conv = (gamma.diff>g.co)
        
        iter.co = sum(iter.co,1)
        if(iter.co>iter.cut){
          warn.clust = 1
          break
        }
      })
  }
  
  #----------------------------------------------------#
  # If optimType = "nlminb" and there is an error, revert
  # to the less sophisticated BFGS and give that a whirl.

  if(class(TroubleInd)=="try-error"&optimType=="nlminb"){
    #--- Initialized ---#
    cdata[["mix"]][["p.est"]] = matrix(test.init,nrow=1,ncol=length(cellTypes))
    cdata[["mix"]][["tau.est"]] = cdata[["beta.est"]][1,]/colSums(cdata[["beta.est"]])
    cdata[["mix"]][["gamma.est"]] = t(t(cdata[["alpha.est"]])/colSums(cdata[["alpha.est"]]))[-c(cdata$I),]/cdata$l_tilde[-cdata$I]
    cdata[["mix"]][["pm.rds.exons"]] = matrix(0,nrow=(cdata$E),ncol=length(cellTypes))
    
    TroubleInd<-try(
      while(any(p.conv,t.conv,g.conv)){
        #--------------------------------------------------#
        # STEP 1 (E-Step): UPDATE the Posterior Means      #
        #--------------------------------------------------#
        cdata[["mix"]][["pm.rds.exons"]] = compute.pmeans(X_fin = cdata$X.fin,
                                                          mix_list = cdata[["mix"]])
        
        #--------------------------------------------------#
        # STEP 2.A (M-Step): UPDATE tau and proportions    #
        #--------------------------------------------------#
        tmp.tp.out = lab_fit_update_v2(cdata = cdata,
                                       nk = length(cellTypes))
        
        #--------------------------------------------------#
        # STEP 2.B (M-Step): UPDATE the Gamma values       #
        #--------------------------------------------------#
        alphas.prev.rds = rbind(cdata[["alpha.est"]],
                                cdata[["mix"]][["gamma.est"]],
                                cdata[["mix"]][["pm.rds.exons"]])
        
        new.gamma = apply(X = alphas.prev.rds,MARGIN = 2,FUN = update.iso.STG2,
                          X.fin = cdata$X.fin, X.prime = cdata$X.prime,
                          I = cdata$I, l_tilde = cdata$l_tilde,optimType = "BFGSWB")
        
        #--------------------------------------------------#
        # CHECK Convergence of the Values / STORE          #
        #--------------------------------------------------#
        #----- Inflation Factors --------------------------#
        # Computing % change causes problems if p and tau have zero entries
        # Best to add inflation factor to control for this problem. 
        
        # P threshold chosen since this gets into the 0.001%
        # In mixture sample. We don't much care about that.
        
        # Tau threshold chosen since at this level we would 
        # average less than 1 read per 10,000,000 total reads
        # in the cluster. Seems a sufficiently small bump. 1 
        # read is approximately zero in my book because of error.
        
        t.inflation = 1e-8
        
        #----- Convergence Computations -------------------#
        new.p = c(tmp.tp.out$p.est,1-sum(tmp.tp.out$p.est))
        p.diff = max(abs(new.p-cdata[["mix"]]$p.est))
        cdata[["mix"]][["p.est"]] = matrix(c(tmp.tp.out[["p.est"]],1-sum(tmp.tp.out[["p.est"]])),nrow=1,ncol=length(cellTypes))
        p.conv = (p.diff>p.co)
        
        tau.diff = max(abs(tmp.tp.out$t.est-cdata[["mix"]]$tau.est))
        cdata[["mix"]][["tau.est"]] = tmp.tp.out[["t.est"]]
        t.conv = (tau.diff>t.co)
        
        gamma.diff = max(abs(cdata[["mix"]][["gamma.est"]]-new.gamma)*cdata$l_tilde[-cdata$I])
        cdata[["mix"]][["gamma.est"]] = new.gamma
        g.conv = (gamma.diff>g.co)
        
        iter.co = sum(iter.co,1)
        if(iter.co>iter.cut){
          warn.clust = 1
          break
        }
      })
  }
  
  #----------------------------------------------------#
  # If optimType = "nlminb" and there is an error, revert
  # to the less sophisticated BFGSWB and give that a whirl.
  
  if(class(TroubleInd)=="try-error"){
    #--- Initialized ---#
    cdata[["mix"]][["p.est"]] = matrix(test.init,nrow=1,ncol=length(cellTypes))
    cdata[["mix"]][["tau.est"]] = cdata[["beta.est"]][1,]/colSums(cdata[["beta.est"]])
    cdata[["mix"]][["gamma.est"]] = t(t(cdata[["alpha.est"]])/colSums(cdata[["alpha.est"]]))[-c(cdata$I),]/cdata$l_tilde[-cdata$I]
    cdata[["mix"]][["pm.rds.exons"]] = matrix(0,nrow=(cdata$E),ncol=length(cellTypes))
    
    TroubleInd<-try(
      while(any(p.conv,t.conv,g.conv)){
        #--------------------------------------------------#
        # STEP 1 (E-Step): UPDATE the Posterior Means      #
        #--------------------------------------------------#
        cdata[["mix"]][["pm.rds.exons"]] = compute.pmeans(X_fin = cdata$X.fin,
                                                          mix_list = cdata[["mix"]])
        
        #--------------------------------------------------#
        # STEP 2.A (M-Step): UPDATE tau and proportions    #
        #--------------------------------------------------#
        tmp.tp.out = lab_fit_update_v2(cdata = cdata,
                                       nk = length(cellTypes))
        
        #--------------------------------------------------#
        # STEP 2.B (M-Step): UPDATE the Gamma values       #
        #--------------------------------------------------#
        alphas.prev.rds = rbind(cdata[["alpha.est"]],
                                cdata[["mix"]][["gamma.est"]],
                                cdata[["mix"]][["pm.rds.exons"]])
        
        new.gamma = apply(X = alphas.prev.rds,MARGIN = 2,FUN = update.iso.STG2,
                          X.fin = cdata$X.fin, X.prime = cdata$X.prime,
                          I = cdata$I, l_tilde = cdata$l_tilde,optimType = "BFGSWB")
        
        #--------------------------------------------------#
        # CHECK Convergence of the Values / STORE          #
        #--------------------------------------------------#
        #----- Inflation Factors --------------------------#
        # Computing % change causes problems if p and tau have zero entries
        # Best to add inflation factor to control for this problem. 
        
        # P threshold chosen since this gets into the 0.001%
        # In mixture sample. We don't much care about that.
        
        # Tau threshold chosen since at this level we would 
        # average less than 1 read per 10,000,000 total reads
        # in the cluster. Seems a sufficiently small bump. 1 
        # read is approximately zero in my book because of error.
        
        t.inflation = 1e-8
        
        #----- Convergence Computations -------------------#
        new.p = c(tmp.tp.out$p.est,1-sum(tmp.tp.out$p.est))
        p.diff = max(abs(new.p-cdata[["mix"]]$p.est))
        cdata[["mix"]][["p.est"]] = matrix(c(tmp.tp.out[["p.est"]],1-sum(tmp.tp.out[["p.est"]])),nrow=1,ncol=length(cellTypes))
        p.conv = (p.diff>p.co)
        
        tau.diff = max(abs(tmp.tp.out$t.est-cdata[["mix"]]$tau.est))
        cdata[["mix"]][["tau.est"]] = tmp.tp.out[["t.est"]]
        t.conv = (tau.diff>t.co)
        
        gamma.diff = max(abs(cdata[["mix"]][["gamma.est"]]-new.gamma)*cdata$l_tilde[-cdata$I])
        cdata[["mix"]][["gamma.est"]] = new.gamma
        g.conv = (gamma.diff>g.co)
        
        iter.co = sum(iter.co,1)
        if(iter.co>iter.cut){
          warn.clust = 1
          break
        }
      })
  }
  
  #----------------------------------------------------#
  # STORE the output                                   #
  #----------------------------------------------------#
  cdata[["CellType_Order"]] = cellTypes
  if(class(TroubleInd)=="try-error"){
    cdata[["WARN"]] = 4
  } else {
    cdata[["WARN"]] = warn.clust
  }
  return(cdata)
}

#------------------------------------------------------#
# WARNING INDICATORS                                   #
#------------------------------------------------------#
# 0 - Optimization Complete
# 1 - Iteration Limit Reached
# 4 - Error in Optimization Routine
# 5 - Optimization not conducted (Error in pure sample fit)

#' @import alabama
STG.Update_Cluster.All<-function(all_data,cellTypes,optimType="nlminb",simple.Init,initPts){
  tmp.out = lapply(X = all_data,FUN = STG.Update_Cluster.SingMI,cellTypes= cellTypes,
                   p.co = 1e-7,t.co=1e-8,g.co=1e-7,optimType=optimType,iter.co=15000,
                   simple.Init=simple.Init, initPts = initPts)
  return(tmp.out)
}


#-------------------------------------------------------------------#
# PROGRAM LIST                                                      #
#-------------------------------------------------------------------#
# (1) STG.Update_Cluster.All -
#     [I.1] all_data    = List item with one element per gene
#               - $CTA       - Pure Data information regarding CTA
#               - $CTB       - Pure Data information regarding CTB
#               - $alpha.est - hyperparams governing gene expression
#               - $beta.est  - hyperparams governing iso expression 
#               - $mix       - Mixture Data information
#                 - $rds_exons_t  - Vector of length E+1 containing reads in mixture
#                 - $gamma.est    - Vector of length I containing iso estimates   
#                 - $tau.est      - parameter of length 1 containing gene estimates
#                 - $p.est        - vector of length K containing estimated proportions
#                 - $pm.rds.exons - Matrix of size ExK containing posterior means
#     [I.2] cellTypes   = vector containing names of cell types used to label
#                         data in the all_data object
#     [I.3] optimType   = nlminb, BFGS, BFGSWB, NM
#     [I.4] simple.Init = TRUE/FALSE whether to use a pre-initialize procedure
#     [I.5] initPts     = a vector (or matrix) of initial points to use for 
#                         the optimization routine.
#
#     Description : Calls to update functions for each element of all_data
#                   Returns list item (1 element per gene), containing relevant
#                   estimates in $mix
#
#########################
# (2) STG.Update_Cluster.SingMI
#     [I.1] cdata       = List item of length 1 containing 1 element of all_data 
#                         ($CTA,$CTB,$mix,pure estimates)
#     [I.2] cellTypes   = vector containing cell type labels used to index all_data
#     [I.3] x.co        = Convergence cutoffs for props, taus (gene), gammas (iso)
#     [I.4] optimType   = nlminb, BFGS, BFGSWB, NM 
#     [I.5] iter.co     = number of iterations allowed
#     [I.6] simple.Init = (TRUE/FALSE) using simple init to decide which initial value
#                         to use
#     [I.7] initPts     = A vector (or matrix) containing a list of start points for K-1 cell types
# 
#     Description : Implements single fit several times across different startpoints to ensure the
#                   best fit possible for our optimization technique.
#
#########################
# (3) applyMI
#       [I.1] x         = a vector of length K containing the start point for the procedure
#       [I.2] cdata     = a single element of the all_data object
#       [I.3] cellTypes = a vector containing the cell type labels used to index cdata
#       [I.4] x.co      = Convergence cutoffs for props, taus (genes), gammas (isos)
#       [I.5] optimType = nlminb, BFGS, BFGSWB, NM
#       [I.6] iter.cut  = Iteration cutoff for the EM algorithm
#
#       Description: This programs applies the initial fit across the grid. This fit
#                    could be simplified so that we estimate an initial likelihood using
#                    initial gene and isoform expressions and certain proportions. 
#
#########################
# (4) STG.Update_Cluster.Single
#       [I.1] cdata     = Single element of all_data to be optimized
#       [I.2] cellTypes = vector containing cell type labels used to index cdata
#       [I.3] x.co      = Convergence criteria for props, taus (gene), and isos (gamma)
#       [I.4] optimType = nlminb, BFGS, BFGSWB, NM
#       [I.5] test.init = a vector of length K containing the starting point of the procedure
#       [I.6] iter.cut  = the maximum number of iterations allowed
#
#       Description: Constructs and outputs results from a single gene fitting. Also provides
#                    a warning indicator with the following output codes
#                    $ WARN
#                       0 - Optimization Complete
#                       1 - Iteration Limit Reached
#                       4 - Error in Optimization Routine
#                       5 - Optimization not conducted (Error in pure sample fit)
#
#########################
# (5) lab_fit_update_v2
#         [I.1] cdata   = data from a single element of all_data
#         [I.2] nk      = a single value specifying the number of cell types
#
#         Description: Performs a single iteration's update of the tau and proportion
#                      estimates returns a list object with $t.est and $p.est
#
#########################
# (6) update.iso.STG2
#          [I.1] X.fin   = Edited design matrix for new gamma parameters (X_i-[l_i/l_I]X_I,...,X_I/l_I)
#          [I.2] X.prime = First (I-1) columns of X.fin pertaining to gamma parameters
#          [I.3] I       = Number of Isoforms used by the considered gene
#          [I.4] l_tilde = Total effective length for an isoform
#          [I.5] alphas.prev.rds = vector of length (2*I-1+E) containing 
#                                  {1:I}           Alpha hyperparameters
#                                  {(I+1):(2*I-1)} previous gamma estimates
#                                  {(2*I):(END)}   posterior mean reads for the optimized cell type
#          [I.6] optimType = nlminb, BFGS, BFGSWB, NM
#   
#          Description: Performs the isoform update procedure for a single cell type, single iteration.
# 
#########################          

