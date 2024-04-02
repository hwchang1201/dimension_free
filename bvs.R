rm(list = ls())

# Rscript bvs.R "moderate" "informed" "yes" 20
set.seed(1234)
args <- commandArgs(TRUE)

cov <- args[1]
method <- args[2]
clip <- args[3]
num_var <- as.integer(args[4]) #how many variables for the intial model
#cov = 'high', 'moderate' 
#method = 'single', 'informed'
#clip = "yes", "no"
#num_var in [0,200]

#cov = "high"
#method = "informed"
#clip = "yes"
#num_var = 20

source("method/rwmh.R")
source("method/imh.R")

# Data generation
n = 200; p = 500; s = 5
b0 = c(8,-12,8,8,-12)
beta = c(b0, rep(0, p - s)) * sqrt( log(p)/n )  # noise sd = 1

C = diag(p)
for (i in 1:p){
  for (j in 1:p){
    if(i != j){
      if (cov == "high"){
        C[i,j] = exp(-abs(i - j)/4)
      }else{
        C[i,j] = exp(-2*abs(i - j))
      }
    }
  }
}
R = chol(C)
# number of iterations
if(method == 'single'){nmcmc=10000}else{nmcmc=1500}

   
X = matrix(rnorm(n*p), nrow=n) %*% R; y = X %*% beta +  rnorm(n)

preprocessed = list(
  XtX = crossprod(X),
  Xty = crossprod(X, y),
  yty = crossprod(y)
)

# true log posterior prob.
g = p^3
kappa = 1
SSRtrue = preprocessed$yty - g/(g+1)*preprocessed$Xty[1:s]%*%solve(preprocessed$XtX[1:s,1:s],preprocessed$Xty[1:s])
truelogpost = lgamma(n/2)- (n/2)*log(pi) - s/2*log(1+g) - n/2*log(SSRtrue) - kappa*s*log(p)
truelogpost = as.numeric(truelogpost)

gammainit = rep(0, p)
gammainit[sample(1:p, num_var )] = 1

if (method == 'single'){

  fit_single = rwmh_bvs(y, X, g, kappa, s0= n,
                             burn = 0, nmc = nmcmc, thin = 1,
                             gammainit = gammainit,
                             verbose = F, debug = F, preprocessed = preprocessed,
                             truegamma = as.logical(beta))
}else{
  if(clip == "yes"){
    fit_single = imh_bvs(y, X, g, kappa, s0= n,
                                          burn = 0, nmc = nmcmc, thin = 1,
                                          gammainit = gammainit,
                                          verbose = F, debug = F, preprocessed = preprocessed,
                                          truegamma = as.logical(beta))
  }else{
    fit_single = imh_bvs_no_clip(y, X, g, kappa, s0= n,
                               burn = 0, nmc = nmcmc, thin = 1,
                               gammainit = gammainit,
                               verbose = F, debug = F, preprocessed = preprocessed,
                               truegamma = as.logical(beta))
  }

  
}


cat(paste("The number of iterations to hit the true: ", fit_single$hit_iter, "/", nmcmc, "\n",
          "The number of iterations to hit the true: ", fit_single$hit_time, "\n",
          "The number of iterations to hit the true: ", fit_single$mcmctime, sep = ""))
