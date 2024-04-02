# THIS CODE IS MADE BY CHANGWOO J. LEE. (c.lee@tamu.edu)
library(mgcv)
#library(microbenchmark)
library(mvnfast) # faster version
library(Matrix)

# auxiliary functions
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


#helper function to compute the marginal probability
SSRgamma_fast <- function(R, g.over.gp1, Xty, yty){
  if(nrow(R)==0){
    Rsq = 0
    SSR = yty # null model
  }else{
    Rsq.numer = sum(forwardsolve(t(R),Xty)^2)
    SSR = yty - g.over.gp1 * Rsq.numer
    Rsq = Rsq.numer / yty
  }
  list(SSR = SSR, Rsq = Rsq)
}



#Yang Wainwright Jordan, 2016, AoS
rwmh_bvs <- function(y, X, g, kappa, s0, burn = 1000,
                          nmc = 5000, thin = 1, gammainit = NULL,
                          verbose = F, debug = F, preprocessed = NULL, 
                          truegamma=NULL){
  # for hitting the true
  found = F
  hit_iter = Inf
  hit_time = Inf
  truegammaidx = which(truegamma==1)
  
  p = ncol(X)
  n = nrow(X)
  gamma = numeric(p)
  if(is.null(gammainit)){
    gamma[sample(1:p, size = floor(s0/2))] <- 1 # random initialize gamma with size s0/2
  }else{
    gamma = gammainit
  }
  gammaidx = which(gamma == 1)
  N=burn+nmc
  effsamp=(N-burn)/thin
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
  if(length(gammaidx)>0){
    Xg = X[,gammaidx, drop = F]
    s = ncol(Xg)
    R = chol(crossprod(Xg))
  }else{
    s = 0
    R = matrix(numeric(0), 0, 0)
  }
  g.over.gp1 = g/(g+1)
  SSRout.current <- SSRgamma_fast(R, g.over.gp1, Xtyoriginal[gammaidx], yty)
  logSSR.current = log(SSRout.current$SSR)
  Rsq = SSRout.current$Rsq
  # MCMC time measure starts 
  mlikelihood_constant = lgamma(n/2) - n*log(pi)/2
  
  t_start = Sys.time()
  record.mcmc = 0
  for(imcmc in 1:N){
    #imcmc = imcmc + 1 
    acc = 0
    #if(runif(1)<0.5){# single flip, do nothing if gamma.prime is out of prior's support 
      idx = sample.int(p, size = 1) # this is scalable(not sensitive with p) 
      if(debug) cat(paste("single filp, s:",s,"\n"))
      # propose Rnew, snew, and gammaidxnew
      if(idx %in% gammaidx){# 1 if delete, 0 if add
        # delete variable
        snew = s - 1
        j = sum(gammaidx <= idx) # can be 1 to s
        Rnew = cholDel2(R, j)# delete one column
        gammaidxnew = gammaidx[-j]
      }else{ # add column
        snew = s + 1
        #if(snew == 1) browser()
        if(snew > s0){#do nothing if snew > s0
          record.mcmc = record.mcmc + 1
          gammaout[[record.mcmc]] <- gammaidx
          accratioout[record.mcmc] <- acc
          #log marginal likelihood
          logpost = mlikelihood_constant - (s/2)*log(1+g)  -s*kappa*log(p) - n/2*logSSR.current
          logpostout[record.mcmc] <- logpost
          next;
        }   
        j = sum(gammaidx <= idx) + 1 # CAN BE 1 TO s+1
        A2 = XtX[idx, idx, drop = F]
        A13 = XtX[gammaidx, idx, drop = F]
        if(j == 1) A1 = NULL else A1 = A13[1:(j-1),,drop = F]
        if(j == s+1) A3 = NULL else A3 = A13[j:s,,drop = F]
        Rnew = cholAdd2(R, j, A2, A1, A3)
        gammaidxnew = append(gammaidx, idx, after = j - 1)
      }
      Xtyselected <- Xtyoriginal[gammaidxnew]
      SSRout.proposal <- SSRgamma_fast(Rnew, g.over.gp1, Xtyselected, yty)
      logSSR.proposal = log(SSRout.proposal$SSR)
      logaccprob = (s - snew)*log(p^kappa*sqrt(1+g)) + (n/2)*(logSSR.current - logSSR.proposal)
      if(log(runif(1))<logaccprob){
        if(debug) cat(paste("single filp accepted, snew:",snew,"\n"))
        acc = 1
        #gamma[idx] <- 1 - gamma[idx]
        gammaidx = gammaidxnew
        s = snew
        R = Rnew
        logSSR.current = logSSR.proposal
        Rsq = SSRout.proposal$Rsq
      }

    if(verbose > 0) if(imcmc %% verbose == 0) {cat(paste("iteration",imcmc,"model size:",s,"\n"))};
    if(imcmc > burn && imcmc%%thin == 0)
    {
      record.mcmc = record.mcmc + 1
      gammaout[[record.mcmc]] <- gammaidx
      accratioout[record.mcmc] <- acc
      #log marginal likelihood
      logpost = mlikelihood_constant - (s/2)*log(1+g) -s*kappa*log(p) - n/2*logSSR.current
      logpostout[record.mcmc] <- logpost
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
  cat("Elapsed time for",N,"MCMC iteration: ",mcmctime,"secs\n")
  
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







onehot <- function(z, n, k){
  Z <- matrix(0,n,k)
  for (j in 1:k){ # not i in 1:n
    Z[which(z==j),j] <- 1
  }
  return(Z)
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



####################################################################################
# GIBBS SAMPLER FOR THE EXTENDED STOCHASTIC BLOCK MODEL  ###########################
####################################################################################

# Inputs:
# Y = VxV symmetric adjacency matrix
# z_init = n-vector of initialization assignment for each node (default = one cluster for each node)
# a,b = parameters of the Beta prior on the thetas

# Output:
# Posterior samples of the community labels for each node i=1,...,n
rwmh_sbm <- function(A, k =5, z_init = NULL, N_iter=10000, a=1, b=1, savez = F, true_z = NULL){
  
  found = F
  hit_iter = Inf
  hit_time = Inf
  if(!is.null(true_z)){
    true_z = cluster_sort(true_z)
    true_nclust = tabulate(true_z)
    true_nclust_sorted = sort(true_nclust)
  }
  
  # singletry
  #ntry= 1
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  n <- nrow(A)
  mode(A) = "logical" # fast indexing
  
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
  #appending_idx1 = rep((0:(ntry-1))*(k*(k+1)/2), each = k)
  #appending_idx2 = (0:(ntry-1))*(k*(k+1)/2)
  
  # ----------------------------------------------
  # Beginning of the M-H sampler
  # ----------------------------------------------
  one_to_k = 1:k
  
  
  t_start = Sys.time()
  for (imcmc in 1:N_iter){
    # ntry = 1
    
    # 1. choose an index j uniformly at random
    j = sample.int(n, size = 1)
    j_oldclust = z[j]
    #browser()
    # cluster number change is not allowed
    if(nclust[j_oldclust]==1){
      print("singleton!")
      next; 
    }
    
    # 2. randomly assign new label to get new assignment
    j_newclust = sample(one_to_k[-j_oldclust], size =1)
    
    
    
    # 3. calculate likelihood ratio
    
    #r_v = crossprod(Z, A[,j]) # same as 
    r_v = tabulate(z[which(A[,j])], nbins = k)
    
    oldidx = index_lookup[,j_oldclust]
    newidx = index_lookup[,j_newclust]
    
    # m_vector
    m_vector_new = m_vector
    m_vector_new[oldidx] = m_vector_new[oldidx] - r_v
    m_vector_new[newidx] = m_vector_new[newidx] + r_v
    
    # N_vector
    N_vector_new = N_vector
    N_vector_new[oldidx] = N_vector_new[oldidx] - nclust
    N_vector_new[newidx] = N_vector_new[newidx] + nclust
    N_vector_new[index_lookup[j_oldclust,j_oldclust]] = N_vector_new[index_lookup[j_oldclust,j_oldclust]] + 1
    N_vector_new[index_lookup[j_oldclust,j_newclust]] = N_vector_new[index_lookup[j_oldclust,j_newclust]] - 1
    
    # comparison..
    # nclust_new = nclust
    # nclust_new[j_oldclust] = nclust_new[j_oldclust]-1
    # nclust_new[j_newclust] = nclust_new[j_newclust]+1
    # 
    # N_full_new = tcrossprod(nclust_new)
    # diag(N_full_new) = nclust_new*(nclust_new-1)/2
    # N_full_new
    
    log_lik_new = sum(lbeta(m_vector_new + a, N_vector_new - m_vector_new + b))
    
    accept = F
    if(log(runif(1)) < log_lik_new - log_lik){
      accept = T
      log_lik = log_lik_new  
      # update z
      
      z[j] <- j_newclust
      #z = GPPM::cluster_sort(z) don't need this
      # update Z (one-hot encoding matrix)
      #Z = onehot(z, n, k)
      
      # update nclust
      #nclust = tabulate(z)}, #maybe simplified
      nclust[j_newclust] = nclust[j_newclust] + 1
      nclust[j_oldclust] = nclust[j_oldclust] - 1
      
      # update m vector, n vector
      m_vector = m_vector_new
      N_vector = N_vector_new
    }
    #browser()
    # store cluster assignments at time imcmc in matrix z_post s.t.
    # z_post[i,t]=h if node i is in cluster h at iteration t
    if(savez) z_post[imcmc,] <- z
    loglik_post[imcmc] <- log_lik
    accept_post[imcmc] <- accept
    #print(table(z_post[,t])) 
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


cluster_sort <- function(x){ 
  nunique = length(unique(x))
  temp = factor(x, labels = 1:nunique)
  temp = as.numeric(as.character(temp))
  res = replace(temp, unique(temp), 1:nunique)[temp]
  res
}


