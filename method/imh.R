# THIS CODE IS MODIFIED VERSION OF THE PREVIOUS CODE BY CHANGWOO J. LEE. (c.lee@tamu.edu)
# data generating function
library(mgcv)
library(microbenchmark)
library(mvnfast) # faster version
library(Matrix)


cholDel2 <- function(R, j){
  n = nrow(R)
  if(j == n){
    R_new = R[-n, -n, drop = F]
  }else{
    R_new = mgcv::choldrop(R, j)
  }
  R_new
}



cholAdd2 <- function(R, j, A2, A1 = NULL, A3 = NULL) {
  n = nrow(R)
  if(j == n+1) {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(n > 0){
      R_new[1:n,1:n] = R
      S12 = drop(backsolve(R, A1, transpose = T))
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[1:(n+1), n+1] = c(S12, S22)
    }else{
      R_new[1,1] = sqrt(as.numeric(A2)) ### fixed - 20210331 / 20210412
    }
  } else {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(j > 1) {
      R11 = R[1:(j-1), 1:(j-1), drop = F]
      R_new[1:(j-1), 1:(j-1)] = R11
      S12 = backsolve(R11, A1, transpose = T)
      R_new[1:(j-1), j] = S12
      S13 = R[1:(j-1), j:n, drop = F]
      R_new[1:(j-1), (j+1):(n+1)] = S13
      
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[j, j] = S22
      S23 = (t(A3) - crossprod(S12, S13)) / S22
      R_new[j, (j+1):(n+1)] = S23
      S33 = mgcv::cholup(R[j:n, j:n, drop = F], S23, up = FALSE) # downdate
      R_new[(j+1):(n+1), (j+1):(n+1)] = S33
    } else {
      S22 = sqrt(as.numeric(A2))
      R_new[1, 1] = S22
      S23 = as.numeric(A3) / S22
      R_new[1, 2:(n+1)] = S23
      S33 = mgcv::cholup(R, S23, up = FALSE) # downdate
      R_new[2:(n+1), 2:(n+1)] = S33
    }
  }
  return(R_new)
}



#Yang Wainwright Jordan, 2016, AoS
# modified to including only the single flip 
# multiple try metropolis - 
# it proposes K proposal index, select 1 with prob. proportional to the posterior 
# then also proposes K-1 auxiliary variables to select it. See Multiple-try Metropolis paper for details
# square root weighting 
imh_bvs <- function(y, X, g, kappa, s0, burn = 1000, nmc = 5000, thin = 1,
                             gammainit = NULL, 
                             verbose = F, debug = F, preprocessed = NULL, K = 20,
                             truegamma = NULL){
  p = ncol(X)
  n = nrow(X)
  K = p
  
  coef_l = p^(1)
  coef_L = p^3
  
  gamma = numeric(p)
  # for hitting the true
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(is.null(gammainit)){
    gamma[sample(1:p, size = floor(s0/2))] <- 1 # random initialize gamma with size s0/2
  }else{
    gamma = gammainit
  }
  gammaidx = which(gamma==1)
  truegammaidx = which(truegamma==1)
  
  nmcmc=burn+nmc
  effsamp=(nmcmc-burn)/thin
  ## output ##
  #betaout=matrix(0,p,effsamp)
  gammaout = list()
  accratioout = numeric(effsamp)
  logpostout = numeric(effsamp)
  Rsqout = numeric(effsamp)
  #sigmaSqout=rep(0,effsamp)
  #precal
  if(is.null(preprocessed)){
    cat("preprocessing XtX, Xty and yty...")
    Xtyoriginal = crossprod(X,y)
    #XtX = as(crossprod(X), "dspMatrix")
    XtX = crossprod(X)
    yty = sum(y^2)
    cat("done\n")
  }else{
    Xtyoriginal = preprocessed$Xty
    #XtX = as(preprocessed$XtX, "dspMatrix")
    XtX = preprocessed$XtX
    yty = as.numeric(preprocessed$yty)
  }
  diagXtX <- diag(XtX)
  if(length(gammaidx)>0){
    Xg = X[,gammaidx, drop = F]
    s = ncol(Xg)
    R = chol(crossprod(Xg))
  }else{
    s = 0
    R = matrix(numeric(0), 0, 0)
  }
  g.over.gp1 = g/(g+1)
  if(nrow(R)==0){
    Rsq = 0
    SSR = yty # null model
  }else{
    bsol = backsolve(R, Xtyoriginal[gammaidx], transpose = T)
    Rsq.numer = sum(bsol^2)
    SSR = yty - g.over.gp1 * Rsq.numer
    Rsq = Rsq.numer / yty
  }
  logSSR.current = log(SSR)
  SSR.current = SSR
  # MCMC time measure starts 
  likelihood_constant = lgamma(n/2) - n*log(pi)/2
  likelihood_constant2 = log(1+g)/2 + kappa*log(p)
  logpost.current = likelihood_constant - s*likelihood_constant2 - (n/2)*log(SSR.current)
  t_start = Sys.time()
  record.mcmc = 0
  for(imcmc in 1:nmcmc){
    #imcmc = imcmc + 1 
    cat(imcmc, "th iteration\n")
    acc = 0
    # Step 1. K independent single flip proposals
    idx_K = 1:p
    #idx_K = sample.int(p, size = K, replace = T) # this is scalable(not sensitive with p)
    #idx_K is not sorted - should we sort?
    if(debug) cat(paste("single filp, s:",s,"\n"))
    
    # Delete variables
    Kidx_drop = which(idx_K%in%gammaidx) # indicies of idx_K that is included in the gammaidx 
    # Add variables
    Kidx_add = which(!(idx_K%in%gammaidx))
    
    logweights = numeric(K)
    logweights_star = numeric(K)
    
    logpost_multiple = numeric(K) # placeholder to save log posterior probabilities for each K
    
    logconstant_drop = likelihood_constant - (s-1)*likelihood_constant2
    
    Rnew_list <- list()
    
    # logpost_multiple for delete. For loop. TODO: can it be done in matrix multiplication?
    for(jj in Kidx_drop){
      j = sum(gammaidx <= idx_K[jj]) # can be 1 to s
      #j = sum(supp <= idx_K[jj]) # can be 1 to s
      Rnew_list[[jj]] <- Rnew <- cholDel2(R, j)# delete one column
      gammaidxnew = gammaidx[-j]
      if(nrow(Rnew)==0){
        SSR = yty # null model
      }else{
        Rsq.numer = sum(backsolve(Rnew, Xtyoriginal[gammaidxnew], transpose = T)^2)
        SSR = yty - g.over.gp1 * Rsq.numer
      }
      logpost_multiple[jj] = logconstant_drop - (n/2)*log(SSR)
      logweights[jj] = (logpost_multiple[jj] - logpost.current)
    }
    
    # logpost_multiple for add.
    logconstant_add = likelihood_constant - (s+1)*likelihood_constant2
    if( s+1 > s0){#do nothing if snew > s0
      logweights[Kidx_add] = -Inf
    }else if(length(Kidx_add)==0){ # no add proposal 
      # do nothing
    }else{
      if(s > 0){
        jvec = idx_K[Kidx_add]
        n.kadd = length(Kidx_add)
        S12.parallel = matrix(backsolve(R, XtX[gammaidx, jvec, drop = F], transpose = T), ncol = n.kadd) # s x n.kadd matrix
        S22.parallel = sqrt(as.numeric(diagXtX[jvec] - colSums(S12.parallel^2))) # length n.kadd vector
        if(any(is.na(S22.parallel))) browser()
        bsol_norm_diff = ((Xtyoriginal[jvec] - crossprod(S12.parallel,bsol))/S22.parallel)^2 # length n.kadd vector
        
        SSR.proposal = SSR.current - g.over.gp1*bsol_norm_diff # amount SSR decreasing, length n.kadd vector
        logpost_multiple[Kidx_add] = logconstant_add - (n/2)*log(SSR.proposal)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
        #browser()
      }else{ # add from null model
        jvec = idx_K[Kidx_add]
        # no S12
        S22.parallel = sqrt(as.numeric(diagXtX[jvec])) # length kadd vector
        bsol_norm_diff = (Xtyoriginal[jvec]/S22.parallel)^2
        SSR.proposal = SSR.current - g.over.gp1*bsol_norm_diff 
        logpost_multiple[Kidx_add] = logconstant_add - (n/2)*log(SSR.proposal)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
      }
    #if(any(is.na(K_weight_normalized))) browser()
      logweights[Kidx_add] = (logpost_multiple[Kidx_add] - logpost.current)
    }
    if(s+1> s0){
      logweights[!Kidx_add] = pmin(pmax(log(coef_l), (logpost_multiple - logpost.current)), log(coef_L))
    }else{
      logweights = pmin(pmax(log(coef_l), (logpost_multiple - logpost.current)), log(coef_L))
    }
    accratio.numer = matrixStats::logSumExp(logweights) #log scale
    weights_normalized = exp(logweights -  accratio.numer)
    
    proposed.Kidx = sample.int(K, size = 1, prob = weights_normalized) # final proposal index of 1:K
    proposed.idx = idx_K[proposed.Kidx] #FINAL PROPOSAL index
    
    #############################################
    log_acc_rate = - log( weights_normalized[proposed.idx])
    
    
    
    
    
    
    
    if(proposed.Kidx %in% Kidx_drop){# drop is accepted
      j = sum(gammaidx <= proposed.idx)
      Rnew = Rnew_list[[proposed.Kidx]]
      gammaidxnew = gammaidx[-j]
      if(length(gammaidxnew)>0){
        bsolnew = backsolve(Rnew,Xtyoriginal[gammaidxnew], transpose = T)
        Rsq.numer = sum(bsolnew^2)
        SSR.proposed = yty - g.over.gp1 * Rsq.numer
        Rsqnew = Rsq.numer / yty
      }else{ # null model
        bsolnew = 0
        SSR.proposed = yty
        Rsqnew = 0
      }
      snew = s-1
    }else{ # add is accepted
      j = sum(gammaidx <= proposed.idx) + 1 # CAN BE 1 TO s+1
      A2 = diagXtX[proposed.idx]
      A13 = XtX[gammaidx, proposed.idx, drop = F]
      if(j == 1) A1 = NULL else A1 = A13[1:(j-1),,drop = F]
      if(j == s+1) A3 = NULL else A3 = A13[j:s,,drop = F]
      Rnew = cholAdd2(R, j, A2, A1, A3)
      gammaidxnew = append(gammaidx, proposed.idx, after = j - 1)
      bsolnew = backsolve(Rnew,Xtyoriginal[gammaidxnew], transpose = T)
      Rsq.numer = sum(bsolnew^2)
      SSR.proposed = yty - g.over.gp1 * Rsq.numer
      Rsqnew = Rsq.numer / yty
      snew = s+1
    }
    
    #min(max(log(coef_l) ,logpost_multiple[proposed.Kidx] - logpost.current), log(coef_L))
    
    logpost.proposed = logpost_multiple[proposed.Kidx]
    #############################################
    
    idx_K_star = 1:p# auxiliary 
    
    Kidx_drop_star = which(idx_K_star%in%gammaidxnew) # indicies of idx_K that is included in the gammaidx 
    Kidx_add_star = which(!(idx_K_star%in%gammaidxnew))
    logpost_multiple_star = numeric(p) # length K -  save log posterior probabilities
    
    logconstant_drop = likelihood_constant - (snew-1)*likelihood_constant2
    Rnew_list <- list()
    
    for(jj in Kidx_drop_star){
      j = sum(gammaidxnew <= idx_K_star[jj]) # can be 1 to s
      Rnew_star <- cholDel2(Rnew, j)# delete one column
      gammaidxnew_star = gammaidxnew[-j]
      if(nrow(Rnew_star)==0){
        SSR_star = yty # null model
      }else{
        Rsq.numer = sum(backsolve(Rnew_star, Xtyoriginal[gammaidxnew_star], transpose = T)^2)
        SSR_star = yty - g.over.gp1 * Rsq.numer
      }
      logpost_multiple_star[jj] = logconstant_drop - (n/2)*log(SSR_star)
      logweights_star[jj] = (logpost_multiple_star[jj] - logpost.proposed)
    }
    
    logconstant_add = likelihood_constant - (snew+1)*likelihood_constant2
    #if( snew+1 > s0){#do nothing if snew > s0
    #  logpost_multiple_star[Kidx_add_star] = -Inf
    #  browser()
      #K_weight_normalized[Kidx_add_star] = 0
    #}else 
    if(length(Kidx_add_star)==0){ # no add proposal 
      # do nothing
    }else{
      if(snew > 0){
        jvec = idx_K_star[Kidx_add_star]
        n.kadd = length(Kidx_add_star)
        S12.parallel = matrix(backsolve(Rnew, XtX[gammaidxnew, jvec, drop = F], transpose = T), ncol = n.kadd) # s x n.kadd matrix
        S22.parallel = sqrt(as.numeric(diagXtX[jvec] - colSums(S12.parallel^2))) # length n.kadd vector
        #if(any(is.na(S22.parallel))) browser()
        if(any(is.na(S22.parallel))) {S22.parallel = rep(1e-5, n.kadd)}
        bsol_norm_diff = ((Xtyoriginal[jvec] - crossprod(S12.parallel,bsolnew))/S22.parallel)^2 # length n.kadd vector
        
        SSR.proposal_star = SSR.proposed - g.over.gp1*bsol_norm_diff # amount SSR decreasing, length n.kadd vector
        logpost_multiple_star[Kidx_add_star] = logconstant_add - (n/2)*log(SSR.proposal_star)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
        #browser()
      }else{ # add from null model
        jvec = idx_K_star[Kidx_add_star]
        # no S12
        S22.parallel = sqrt(as.numeric(diagXtX[jvec])) # length kadd vector
        bsol_norm_diff = (Xtyoriginal[jvec]/S22.parallel)^2
        SSR.proposal_star = SSR.proposed - g.over.gp1*bsol_norm_diff 
        logpost_multiple_star[Kidx_add_star] = logconstant_add - (n/2)*log(SSR.proposal_star)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
      }
      #if(any(is.na(K_weight_normalized))) browser()
      logweights_star[Kidx_add_star] = (logpost_multiple_star[Kidx_add_star] - logpost.proposed)
    }
    
    
    logweights_star = pmin(pmax(log(coef_l), (logpost_multiple_star - logpost.proposed)), log(coef_L))
    log_acc_rate = log_acc_rate  + min(max(log(coef_l) ,logpost.current-logpost_multiple[proposed.Kidx]), log(coef_L)) - matrixStats::logSumExp(logweights_star) + 
      (logpost.proposed - logpost.current)
    
    #########################################################################################
    if(any(is.na(logpost_multiple_star))) browser()
    
    #accratio.denom = matrixStats::logSumExp(c(logweights_star, (logpost.current - logpost.proposed)/2 )) # log scale
    
    if(log(runif(1)) < log_acc_rate){
      acc = 1
      if(debug) cat(paste("single filp accepted,\n"))
      R = Rnew
      gammaidx = gammaidxnew
      bsol = bsolnew
      SSR.current = SSR.proposed
      Rsq = Rsqnew
      s = snew
      logpost.current = logpost.proposed
    }
    
    
    #if(length(gammaidx)==9) return(imcmc)
    #step 1b(optional)draw phi(precision) and step 2(optional)draw beta: omitted
    if(verbose > 0) if(imcmc %% verbose == 0) {cat(paste("iteration",imcmc,"model size:",s,"\n"))};
    if(imcmc > burn && imcmc%%thin == 0)
    {
      record.mcmc = record.mcmc + 1
      gammaout[[record.mcmc]] <- gammaidx
      accratioout[record.mcmc] <- acc
      logpostout[record.mcmc] <- logpost.current
      Rsqout[record.mcmc] <- Rsq
      
      # record hit-time
      if(!is.null(truegamma) & !found){
        if(length(gammaidx)==length(truegammaidx)){
          if(all(sort(gammaidx) == truegammaidx)){
            hit_iter = imcmc
            found = T
            # summarize
            hit_time = difftime(Sys.time(), t_start, units = "secs")
          }
        }
      }
      
    }
  }
  # summarize
  mcmctime = difftime(Sys.time(), t_start, units = "secs")
  cat("Elapsed time for",nmcmc,"MCMC iteration: ",mcmctime,"secs\n")
  
  bin <- tabulate(unlist(gammaout), nbins = p)
  PIP = bin/record.mcmc
  
  list(pip = PIP,
       gammaout = gammaout,
       accratio = mean(accratioout),
       logpostout = logpostout,
       Rsqout = Rsqout, 
       mcmctime = mcmctime,
       hit_time = hit_time,
       hit_iter = hit_iter)
}


imh_bvs_no_clip <- function(y, X, g, kappa, s0, burn = 1000, nmc = 5000, thin = 1,
                          gammainit = NULL, 
                          verbose = F, debug = F, preprocessed = NULL, K = 20,
                          truegamma = NULL){
  p = ncol(X)
  n = nrow(X)
  K = p
  
  coef_l = p^{-200}
  coef_L = p^{200}
  
  gamma = numeric(p)
  # for hitting the true
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(is.null(gammainit)){
    gamma[sample(1:p, size = floor(s0/2))] <- 1 # random initialize gamma with size s0/2
  }else{
    gamma = gammainit
  }
  gammaidx = which(gamma==1)
  truegammaidx = which(truegamma==1)
  
  nmcmc=burn+nmc
  effsamp=(nmcmc-burn)/thin
  ## output ##
  #betaout=matrix(0,p,effsamp)
  gammaout = list()
  accratioout = numeric(effsamp)
  logpostout = numeric(effsamp)
  Rsqout = numeric(effsamp)
  #sigmaSqout=rep(0,effsamp)
  #precal
  if(is.null(preprocessed)){
    cat("preprocessing XtX, Xty and yty...")
    Xtyoriginal = crossprod(X,y)
    #XtX = as(crossprod(X), "dspMatrix")
    XtX = crossprod(X)
    yty = sum(y^2)
    cat("done\n")
  }else{
    Xtyoriginal = preprocessed$Xty
    #XtX = as(preprocessed$XtX, "dspMatrix")
    XtX = preprocessed$XtX
    yty = as.numeric(preprocessed$yty)
  }
  diagXtX <- diag(XtX)
  if(length(gammaidx)>0){
    Xg = X[,gammaidx, drop = F]
    s = ncol(Xg)
    R = chol(crossprod(Xg))
  }else{
    s = 0
    R = matrix(numeric(0), 0, 0)
  }
  g.over.gp1 = g/(g+1)
  if(nrow(R)==0){
    Rsq = 0
    SSR = yty # null model
  }else{
    bsol = backsolve(R, Xtyoriginal[gammaidx], transpose = T)
    Rsq.numer = sum(bsol^2)
    SSR = yty - g.over.gp1 * Rsq.numer
    Rsq = Rsq.numer / yty
  }
  logSSR.current = log(SSR)
  SSR.current = SSR
  # MCMC time measure starts 
  likelihood_constant = lgamma(n/2) - n*log(pi)/2
  likelihood_constant2 = log(1+g)/2 + kappa*log(p)
  logpost.current = likelihood_constant - s*likelihood_constant2 - (n/2)*log(SSR.current)
  t_start = Sys.time()
  record.mcmc = 0
  for(imcmc in 1:nmcmc){
    #imcmc = imcmc + 1 
    cat(imcmc, "th iteration\n")
    acc = 0
    # Step 1. K independent single flip proposals
    idx_K = 1:p
    #idx_K = sample.int(p, size = K, replace = T) # this is scalable(not sensitive with p)
    #idx_K is not sorted - should we sort?
    if(debug) cat(paste("single filp, s:",s,"\n"))
    
    # Delete variables
    Kidx_drop = which(idx_K%in%gammaidx) # indicies of idx_K that is included in the gammaidx 
    # Add variables
    Kidx_add = which(!(idx_K%in%gammaidx))
    
    logweights = numeric(K)
    logweights_star = numeric(p)
    
    logpost_multiple = numeric(K) # placeholder to save log posterior probabilities for each K
    
    logconstant_drop = likelihood_constant - (s-1)*likelihood_constant2
    
    Rnew_list <- list()
    
    # logpost_multiple for delete. For loop. TODO: can it be done in matrix multiplication?
    for(jj in Kidx_drop){
      j = sum(gammaidx <= idx_K[jj]) # can be 1 to s
      Rnew_list[[jj]] <- Rnew <- cholDel2(R, j)# delete one column
      gammaidxnew = gammaidx[-j]
      if(nrow(Rnew)==0){
        SSR = yty # null model
      }else{
        Rsq.numer = sum(backsolve(Rnew, Xtyoriginal[gammaidxnew], transpose = T)^2)
        SSR = yty - g.over.gp1 * Rsq.numer
      }
      logpost_multiple[jj] = logconstant_drop - (n/2)*log(SSR)
      logweights[jj] = (logpost_multiple[jj] - logpost.current)
    }
    
    # logpost_multiple for add.
    logconstant_add = likelihood_constant - (s+1)*likelihood_constant2
    if( s+1 > s0){#do nothing if snew > s0
      logweights[Kidx_add] = -Inf
    }else if(length(Kidx_add)==0){ # no add proposal 
      # do nothing
    }else{
      if(s > 0){
        jvec = idx_K[Kidx_add]
        n.kadd = length(Kidx_add)
        S12.parallel = matrix(backsolve(R, XtX[gammaidx, jvec, drop = F], transpose = T), ncol = n.kadd) # s x n.kadd matrix
        S22.parallel = sqrt(as.numeric(diagXtX[jvec] - colSums(S12.parallel^2))) # length n.kadd vector
        if(any(is.na(S22.parallel))) browser()
        bsol_norm_diff = ((Xtyoriginal[jvec] - crossprod(S12.parallel,bsol))/S22.parallel)^2 # length n.kadd vector
        
        SSR.proposal = SSR.current - g.over.gp1*bsol_norm_diff # amount SSR decreasing, length n.kadd vector
        logpost_multiple[Kidx_add] = logconstant_add - (n/2)*log(SSR.proposal)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
        #browser()
      }else{ # add from null model
        jvec = idx_K[Kidx_add]
        # no S12
        S22.parallel = sqrt(as.numeric(diagXtX[jvec])) # length kadd vector
        bsol_norm_diff = (Xtyoriginal[jvec]/S22.parallel)^2
        SSR.proposal = SSR.current - g.over.gp1*bsol_norm_diff 
        logpost_multiple[Kidx_add] = logconstant_add - (n/2)*log(SSR.proposal)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
      }
      #if(any(is.na(K_weight_normalized))) browser()
      logweights[Kidx_add] = (logpost_multiple[Kidx_add] - logpost.current)
    }
    
    logweights = pmin(pmax(log(coef_l), (logpost_multiple - logpost.current)), log(coef_L))
    accratio.numer = matrixStats::logSumExp(logweights) #log scale
    weights_normalized = exp(logweights -  accratio.numer)
    
    proposed.Kidx = sample.int(K, size = 1, prob = weights_normalized) # final proposal index of 1:K
    proposed.idx = idx_K[proposed.Kidx] #FINAL PROPOSAL index
    
    #############################################
    log_acc_rate = - log( weights_normalized[proposed.idx])
    
    
    
    
    
    
    
    if(proposed.Kidx %in% Kidx_drop){# drop is accepted
      j = sum(gammaidx <= proposed.idx)
      Rnew = Rnew_list[[proposed.Kidx]]
      gammaidxnew = gammaidx[-j]
      if(length(gammaidxnew)>0){
        bsolnew = backsolve(Rnew,Xtyoriginal[gammaidxnew], transpose = T)
        Rsq.numer = sum(bsolnew^2)
        SSR.proposed = yty - g.over.gp1 * Rsq.numer
        Rsqnew = Rsq.numer / yty
      }else{ # null model
        bsolnew = 0
        SSR.proposed = yty
        Rsqnew = 0
      }
      snew = s-1
    }else{ # add is accepted
      j = sum(gammaidx <= proposed.idx) + 1 # CAN BE 1 TO s+1
      A2 = diagXtX[proposed.idx]
      A13 = XtX[gammaidx, proposed.idx, drop = F]
      if(j == 1) A1 = NULL else A1 = A13[1:(j-1),,drop = F]
      if(j == s+1) A3 = NULL else A3 = A13[j:s,,drop = F]
      Rnew = cholAdd2(R, j, A2, A1, A3)
      gammaidxnew = append(gammaidx, proposed.idx, after = j - 1)
      bsolnew = backsolve(Rnew,Xtyoriginal[gammaidxnew], transpose = T)
      Rsq.numer = sum(bsolnew^2)
      SSR.proposed = yty - g.over.gp1 * Rsq.numer
      Rsqnew = Rsq.numer / yty
      snew = s+1
    }
    
    #min(max(log(coef_l) ,logpost_multiple[proposed.Kidx] - logpost.current), log(coef_L))
    
    logpost.proposed = logpost_multiple[proposed.Kidx]
    #############################################
    
    idx_K_star = 1:p# auxiliary 
    
    Kidx_drop_star = which(idx_K_star%in%gammaidxnew) # indicies of idx_K that is included in the gammaidx 
    Kidx_add_star = which(!(idx_K_star%in%gammaidxnew))
    logpost_multiple_star = numeric(p) # length K -  save log posterior probabilities
    
    logconstant_drop = likelihood_constant - (snew-1)*likelihood_constant2
    Rnew_list <- list()
    
    for(jj in Kidx_drop_star){
      j = sum(gammaidxnew <= idx_K_star[jj]) # can be 1 to s
      Rnew_star <- cholDel2(Rnew, j)# delete one column
      gammaidxnew_star = gammaidxnew[-j]
      if(nrow(Rnew_star)==0){
        SSR_star = yty # null model
      }else{
        Rsq.numer = sum(backsolve(Rnew_star, Xtyoriginal[gammaidxnew_star], transpose = T)^2)
        SSR_star = yty - g.over.gp1 * Rsq.numer
      }
      logpost_multiple_star[jj] = logconstant_drop - (n/2)*log(SSR_star)
      logweights_star[jj] = (logpost_multiple_star[jj] - logpost.proposed)
    }
    
    logconstant_add = likelihood_constant - (snew+1)*likelihood_constant2
    #if( snew+1 > s0){#do nothing if snew > s0
    #  logpost_multiple_star[Kidx_add_star] = -Inf
    #  browser()
    #K_weight_normalized[Kidx_add_star] = 0
    #}else 
    if(length(Kidx_add_star)==0){ # no add proposal 
      # do nothing
    }else{
      if(snew > 0){
        jvec = idx_K_star[Kidx_add_star]
        n.kadd = length(Kidx_add_star)
        S12.parallel = matrix(backsolve(Rnew, XtX[gammaidxnew, jvec, drop = F], transpose = T), ncol = n.kadd) # s x n.kadd matrix
        S22.parallel = sqrt(as.numeric(diagXtX[jvec] - colSums(S12.parallel^2))) # length n.kadd vector
        if(any(is.na(S22.parallel))) browser()
        bsol_norm_diff = ((Xtyoriginal[jvec] - crossprod(S12.parallel,bsolnew))/S22.parallel)^2 # length n.kadd vector
        
        SSR.proposal_star = SSR.proposed - g.over.gp1*bsol_norm_diff # amount SSR decreasing, length n.kadd vector
        logpost_multiple_star[Kidx_add_star] = logconstant_add - (n/2)*log(SSR.proposal_star)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
        #browser()
      }else{ # add from null model
        jvec = idx_K_star[Kidx_add_star]
        # no S12
        S22.parallel = sqrt(as.numeric(diagXtX[jvec])) # length kadd vector
        bsol_norm_diff = (Xtyoriginal[jvec]/S22.parallel)^2
        SSR.proposal_star = SSR.proposed - g.over.gp1*bsol_norm_diff 
        logpost_multiple_star[Kidx_add_star] = logconstant_add - (n/2)*log(SSR.proposal_star)
        #K_weight_normalized = exp(logpost_multiple - matrixStats::logSumExp(logpost_multiple))
      }
      #if(any(is.na(K_weight_normalized))) browser()
      logweights_star[Kidx_add_star] = (logpost_multiple_star[Kidx_add_star] - logpost.proposed)
    }
    
    
    logweights_star = pmin(pmax(log(coef_l), (logpost_multiple_star - logpost.proposed)), log(coef_L))
    log_acc_rate = log_acc_rate  + min(max(log(coef_l) ,logpost.current-logpost_multiple[proposed.Kidx]), log(coef_L)) - matrixStats::logSumExp(logweights_star) + 
      (logpost.proposed - logpost.current)
    
    #########################################################################################
    if(any(is.na(logpost_multiple_star))) browser()
    
    #accratio.denom = matrixStats::logSumExp(c(logweights_star, (logpost.current - logpost.proposed)/2 )) # log scale
    
    if(log(runif(1)) < log_acc_rate){
      acc = 1
      if(debug) cat(paste("single filp accepted,\n"))
      R = Rnew
      gammaidx = gammaidxnew
      bsol = bsolnew
      SSR.current = SSR.proposed
      Rsq = Rsqnew
      s = snew
      logpost.current = logpost.proposed
    }
    
    
    #if(length(gammaidx)==9) return(imcmc)
    #step 1b(optional)draw phi(precision) and step 2(optional)draw beta: omitted
    if(verbose > 0) if(imcmc %% verbose == 0) {cat(paste("iteration",imcmc,"model size:",s,"\n"))};
    if(imcmc > burn && imcmc%%thin == 0)
    {
      record.mcmc = record.mcmc + 1
      gammaout[[record.mcmc]] <- gammaidx
      accratioout[record.mcmc] <- acc
      logpostout[record.mcmc] <- logpost.current
      Rsqout[record.mcmc] <- Rsq
      
      # record hit-time
      if(!is.null(truegamma) & !found){
        if(length(gammaidx)==length(truegammaidx)){
          if(all(sort(gammaidx) == truegammaidx)){
            hit_iter = imcmc
            found = T
            # summarize
            hit_time = difftime(Sys.time(), t_start, units = "secs")
          }
        }
      }
      
    }
  }
  # summarize
  mcmctime = difftime(Sys.time(), t_start, units = "secs")
  cat("Elapsed time for",nmcmc,"MCMC iteration: ",mcmctime,"secs\n")
  
  bin <- tabulate(unlist(gammaout), nbins = p)
  PIP = bin/record.mcmc
  
  list(pip = PIP,
       gammaout = gammaout,
       accratio = mean(accratioout),
       logpostout = logpostout,
       Rsqout = Rsqout, 
       mcmctime = mcmctime,
       hit_time = hit_time,
       hit_iter = hit_iter)
}


myf = function(x, z, k){
  tabulate(z[x], nbins = k)
}

imh_sbm<- function(A, k =5, ntry = 10, z_init = NULL, N_iter=10000, a=1, b=1, 
                             savez = F, true_z = NULL, sparse = F){
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  # for hitting the true
  n = nrow(A)
  coef_l = n
  coef_L = n^3
  
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(!is.null(true_z)){
    true_z = cluster_sort(true_z)
    true_nclust = tabulate(true_z)
    true_nclust_sorted = sort(true_nclust)
  }
  
  
  n <- nrow(A)
  if(!sparse){
    #mode(A) = "logical" # fast indexing
  }else{
    A = Matrix::Matrix(A)
    A = as(A, "nMatrix") ## important, improves speed almost x2
  }
  
  Adj_Edgelist = list()
  for(i in 1:n) Adj_Edgelist[[i]] = which(A[,i]==1)
  
  if(is.null(z_init)){
    z_init = rep(1:k, each = n/k, length = n)
    #z_init = sample(z_init)
    #z_init = rep(1:k, times = c(5, 16, 19, 26, 34))
  }else{
    if(length(z_init)!=n) stop("wrong z_init length")
    if(length(unique(z_init))!=k) stop("wrong ncluster of z_init") 
  }
  #z_init = GPPM::cluster_sort(z_init)
  z = z_init
  mode(z) = "integer"
  # cluster assignments are encoded in two equivalent waAs:
  # [i] a VxH matrix Z, s.t. Z[i,h]=1{node i in cluster h}, faster to use within each iteration
  Z <- onehot(z, n, k)
  
  # [ii] a vector of length n containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[i,t]=h if node i is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  if(savez) z_post <- matrix(NA,N_iter, n)
  loglik_post <- numeric(N_iter)
  accept_post <- numeric(N_iter)
  # Create the matrix with block connections
  temp   <- A%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # number of edges between clusters, represented by vector
  m_vector = m_full[upper.tri(m_full, diag = T)]
  
  nclust = tabulate(z)
  
  N_full = tcrossprod(nclust)
  diag(N_full) = nclust*(nclust-1)/2
  # number of possible edges between clusters, represented by vector
  N_vector = N_full[upper.tri(N_full, diag = T)]
  
  log_lik = sum(lbeta(m_vector + a, N_vector - m_vector + b))
  
  
  # Index lookuptable
  index_lookup = matrix(0L, k, k)
  index_lookup[upper.tri(index_lookup, diag = T)] <- 1:(k*(k+1)/2)
  index_lookup = (index_lookup + t(index_lookup))
  diag(index_lookup) = diag(index_lookup)/2
  #index_lookup
  appending_idx1 = rep((0:(n-1))*(k*(k+1)/2), each = k)
  appending_idx2 = (0:(n-1))*(k*(k+1)/2)
  
  appending_idx3 = rep((0:(n-1))*(k*(k+1)/2), each = k)
  appending_idx4 = (0:(n-1))*(k*(k+1)/2)
  
  # ----------------------------------------------
  # Beginning of the M-H sampler
  # ----------------------------------------------
  one_to_k = 1:k
  
  t_start = Sys.time()
  for (imcmc in 1:N_iter){
    
    # j is now vector 
    # j = sample.int(n, size = ntry, replace = T)
    j = 1:n
    j_oldclust = z[j]
    # cluster number change is not allowed
    if(any(nclust[j_oldclust]==1)) next;
    
    # 2-multi. 
    rand_km1 = sample.int(k-1, size = ntry, replace = T)
    rand_km1 = rep(1, n)
    j_newclust = rand_km1 + (j_oldclust<=rand_km1)
    
    # 3. calculate likelihood ratio
    #browser()
    if(sparse){
      r_v = as.vector(crossprod(Z, A[,j])) # same as tabulate(z[which(A[,j])], nbins = k)
    }else{
      r_v = as.vector(vapply(Adj_Edgelist[j], FUN = myf, z = z, k = k, FUN.VALUE = integer(k)))
    }
    # microbenchmark(
    # #r_v = as.vector(crossprod(Z, A[,j])), # same as tabulate(z[which(A[,j])], nbins = k)
    # #r_v2 = as.vector(sapply(Adj_Edgelist[j], FUN = f, z = z)),
    # r_v3 = as.vector(vapply(Adj_Edgelist[j], FUN = myf, z = z, k=k,FUN.VALUE = integer(k))),
    # r_v4 = as.vector(vapply(Adj_Edgelist[j], FUN = myf, z = zint, k=k,FUN.VALUE = integer(k)))
    # )
    
    oldidx = as.vector(index_lookup[,j_oldclust]) + appending_idx1 # appended, vector
    newidx = as.vector(index_lookup[,j_newclust]) + appending_idx1 # appended, vector
    idx_plusone = index_lookup[cbind(j_oldclust,j_oldclust)] + appending_idx2
    idx_minusone = index_lookup[cbind(j_oldclust,j_newclust)]+ appending_idx2
    #browser()
    # each column is each try
    m_vector_appended_new = rep(m_vector, n)
    m_vector_appended_new[oldidx] = m_vector_appended_new[oldidx] - r_v
    m_vector_appended_new[newidx] = m_vector_appended_new[newidx] + r_v
    
    m_matrix_new = matrix(m_vector_appended_new, ncol = n)
    
    # N_vector
    N_vector_appended_new = rep(N_vector, n)
    N_vector_appended_new[oldidx] = N_vector_appended_new[oldidx] - nclust
    N_vector_appended_new[newidx] = N_vector_appended_new[newidx] + nclust
    N_vector_appended_new[idx_plusone] = N_vector_appended_new[idx_plusone] +1
    N_vector_appended_new[idx_minusone] = N_vector_appended_new[idx_minusone] -1
    
    N_matrix_new = matrix(N_vector_appended_new, ncol = n)
    
    # ntry vector
    log_lik_new = colSums(lbeta(m_matrix_new + a, N_matrix_new - m_matrix_new + b))
    
    logweights = pmin(pmax(log(coef_l), (log_lik_new - log_lik)), log(coef_L))
    
    
    # logweights = (log_lik_new - log_lik)/2
    accratio.numer = matrixStats::logSumExp(logweights) #log scale
    weights_normalized = exp(logweights - accratio.numer)
    
    proposed.Kidx = sample(1:n, size = 1, prob = weights_normalized) # final proposal index of 1:K
    proposed.idx = j[proposed.Kidx] # proposed single flip
    log_acc_rate = - log( weights_normalized[proposed.idx])
    
    # next 
    #browser()
    proposed.z = z
    proposed.z[proposed.idx] = j_newclust[proposed.Kidx]
    proposed.m_vector = m_matrix_new[,proposed.Kidx]
    proposed.N_vector = N_matrix_new[,proposed.Kidx]
    proposed.loglik = log_lik_new[proposed.Kidx]
    
    proposed.nclust = nclust
    proposed.nclust[j_newclust[proposed.Kidx]] = proposed.nclust[j_newclust[proposed.Kidx]] + 1
    proposed.nclust[j_oldclust[proposed.Kidx]] = proposed.nclust[j_oldclust[proposed.Kidx]] - 1
    
    if(sparse){
      proposed.Z = Z
      proposed.Z[j[proposed.Kidx], j_oldclust[proposed.Kidx]] = 0
      proposed.Z[j[proposed.Kidx], j_newclust[proposed.Kidx]] = 1
    }
    
    ## single flip from the proposed state
    #jstar = sample.int(n, size = ntry -1 , replace = T)
    jstar = 1:n
    jstar_oldclust = proposed.z[jstar]
    
    # cluster number change is not allowed
    if(any(proposed.nclust[jstar_oldclust]==1)) next;
    
    # 2-multi. 
    #rand_km1 = sample.int(k-1, size = ntry-1, replace = T)
    rand_km1 = rep(1, n)
    jstar_newclust = rand_km1 + (jstar_oldclust<=rand_km1)
    
    
    # 3. calculate likelihood ratio
    if(sparse){
      r_v = as.vector(crossprod(proposed.Z, A[,jstar]))
    }else{
      r_v = as.vector(vapply(Adj_Edgelist[jstar], FUN = myf, z = proposed.z, k = k, FUN.VALUE = integer(k)))
    }
    
    oldidx = as.vector(index_lookup[,jstar_oldclust]) + appending_idx3 # appended, vector
    newidx = as.vector(index_lookup[,jstar_newclust]) + appending_idx3 # appended, vector
    idx_plusone = index_lookup[cbind(jstar_oldclust,jstar_oldclust)] + appending_idx4
    idx_minusone = index_lookup[cbind(jstar_oldclust,jstar_newclust)]+ appending_idx4
    
    # each column is each try
    m_vector_appended_new = rep(proposed.m_vector, n)
    m_vector_appended_new[oldidx] = m_vector_appended_new[oldidx] - r_v
    m_vector_appended_new[newidx] = m_vector_appended_new[newidx] + r_v
    
    m_matrix_newstar = matrix(m_vector_appended_new, ncol = n)
    
    # N_vector
    N_vector_appended_new = rep(proposed.N_vector, n)
    N_vector_appended_new[oldidx] = N_vector_appended_new[oldidx] - proposed.nclust
    N_vector_appended_new[newidx] = N_vector_appended_new[newidx] + proposed.nclust
    N_vector_appended_new[idx_plusone] = N_vector_appended_new[idx_plusone] +1
    N_vector_appended_new[idx_minusone] = N_vector_appended_new[idx_minusone] -1
    
    N_matrix_newstar = matrix(N_vector_appended_new, ncol =n)
    
    # ntry -1 vector
    log_lik_newstar = colSums(lbeta(m_matrix_newstar + a, N_matrix_newstar - m_matrix_newstar + b))
    
    logweights_star = pmin(pmax(log(coef_l), (log_lik_newstar - proposed.loglik)), log(coef_L))
    log_acc_rate = log_acc_rate  + min(max(log(coef_l) ,log_lik - log_lik_new[proposed.idx]) , log(coef_L)) - matrixStats::logSumExp(logweights_star) + 
      (log_lik_new[proposed.idx] - log_lik)
    
    (log_lik_new - log_lik)[proposed.idx]
    accratio.denom = matrixStats::logSumExp(c(logweights_star, -logweights[proposed.Kidx] )) # log scale
    
    accept = F
    #if(log(runif(1)) < accratio.numer - accratio.denom){
    if(log(runif(1)) < exp(log_acc_rate)){
      accept = T
      z = proposed.z
      if(sparse) Z = proposed.Z
      m_vector = proposed.m_vector 
      N_vector = proposed.N_vector
      log_lik = proposed.loglik
      nclust = proposed.nclust 
      
    }
    
    
    if(savez) z_post[imcmc,] <- z
    loglik_post[imcmc] <- log_lik
    accept_post[imcmc] <- accept
    
    # record hit-time
    if(!is.null(true_z) & !found){
      if(all(sort(nclust)==true_nclust_sorted)){ 
        if(all(cluster_sort(z) == true_z)){
          hit_iter = imcmc
          found = T
          # summarize
          hit_time = difftime(Sys.time(), t_start, units = "secs")
        }
      }
    }
    
    #print(table(z_post[,t])) 
    #if (imcmc%%1000 == 0){print(paste("Iteration:", imcmc))}
  }
  mcmctime = difftime(Sys.time(), t_start, units = "secs")
  cat("Elapsed time for",N_iter,"MCMC iteration: ",mcmctime,"secs\n")
  
  out = list()
  if(savez) out$z_post = z_post
  out$loglik_post = loglik_post
  out$accept_post = accept_post
  out$mcmctime = mcmctime
  out$hit_time = hit_time
  out$hit_iter = hit_iter
  return(out)
}





onehot <- function(z, n, k){
  Z <- matrix(0,n,k)
  for (j in 1:k){ # not i in 1:n
    Z[which(z==j),j] <- 1
  }
  return(Z)
}

cluster_sort <- function(x){ 
  nunique = length(unique(x))
  temp = factor(x, labels = 1:nunique)
  temp = as.numeric(as.character(temp))
  res = replace(temp, unique(temp), 1:nunique)[temp]
  res
}

calculatelog_lik <- function(A, z){
  n = length(z)
  k = length(unique(z))
  Z = onehot(z, n, k)
  # Create the matrix with block connections
  mode(A) = "logical"
  temp   <- A%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),k)
  
  # number of edges between clusters, represented by vector
  m_vector = m_full[upper.tri(m_full, diag = T)]
  
  nclust = tabulate(z)
  
  N_full = tcrossprod(nclust)
  diag(N_full) = nclust*(nclust-1)/2
  # number of possible edges between clusters, represented by vector
  N_vector = N_full[upper.tri(N_full, diag = T)]
  
  log_lik = sum(lbeta(m_vector + 1, N_vector - m_vector + 1))
  log_lik
}


