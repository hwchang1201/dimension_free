rm(list = ls())
library(igraph)

set.seed(1234)
args <- commandArgs(TRUE)
method <- args[1]
initial <- args[2]

# Rscript sbm.R "single" "good"

# method = "single", "informed"
# initial = "good", "bad"

source("method/rwmh.R")
source("method/imh.R")

n= 1000; k= 2
p = 10^(-1)   # within-prob
q = 10^(-8) # cross-prob
n*(sqrt(p)-sqrt(q))^2/(k*log(n))  # around 4
# 
# ################################################################



Q = matrix(q, k, k)
diag(Q) = p
g = sample_sbm(n, pref.matrix = Q, block.sizes = c(n/2,n/2))
# g = sample_sbm(n, pref.matrix = Q, block.sizes = c(270,460))
components(g)$no 

true_z = c(rep(1, n/2), rep(2, n/2))
# true_z = c(rep(1, 270), rep(2, 460))

A = as_adj(g, sparse = F)


truelogpost = calculatelog_lik(A, true_z)
if(method == "single"){
  N_iter = 20000
} else{
  N_iter = 2000
}


if(initial=="bad"){
  z_init = true_z
  wrongidx = sample.int(n, size = n/2)
  z_init[wrongidx] = (true_z[wrongidx] + sample(k-1, size = n/2, replace = T))%% k
  z_init[which(z_init==0)] = k
}else{
  z_init = true_z
  wrongidx = sample.int(n, size = n/4)
  z_init[wrongidx] = (true_z[wrongidx] + sample(k-1, size = n/4, replace = T))%% k
  z_init[which(z_init==0)] = k
}



if(method == "single"){
  fit_single = rwmh_sbm(A, k, z_init = z_init,
                             N_iter = N_iter, true_z = true_z, savez = T)
}else{
  fit_single = imh_sbm(A, k, ntry = n, z_init = z_init, 
                                 N_iter = N_iter, true_z = true_z, savez = T)
}


cat(paste("The number of iterations to hit the true: ", fit_single$hit_iter, "/", N_iter, "\n",
          "The number of iterations to hit the true: ", fit_single$hit_time, "\n",
          "The number of iterations to hit the true: ", fit_single$mcmctime, sep = ""))
