rm(list = ls())

setwd("/Users/dehuabi/Desktop/PAMHC_sensitivity/")

library(Rlab)
library(boot)

m <- 1

# Treatment effect
gen_out <- function(n3, n1_actual, n2_actual, delta, cov_d, T_d){
  # y_d_cont, y_d_bin the continuous/binary outcome in G1 and G2
  # y_d_cont_rwd, y_d_bin_rwd the continuous/binary outcomes in G3
  y_d_cont <- rep(NA, (n1_actual+n2_actual))
  for(i in 1:(n1_actual+n2_actual)){
    if(T_d[i] == 1){
      # Treatment arm
      y_d_cont[i] <- sum(cov_d[i,c(1,2)]) + delta + rnorm(1)
    }else{
      # Control arm
      y_d_cont[i] <- sum(cov_d[i,c(1,2)]) + rnorm(1)
    }
  }
  y_d_cont_rwd <- rep(NA, n3)
  for(i in 1:n3){
    # Generate continuous and binary outcomes
    y_d_cont_rwd[i] <- sum(cov_d_rwd[i,c(1,3)]) + rnorm(1)
  }
  
  y_cont_mat <- matrix(rep(NA, max(c(n1_actual, n2_actual, n3))*3), nrow = 3)
  y_cont_mat[1,1:n1_actual] <- y_d_cont[T_d == 1]
  y_cont_mat[2,1:n2_actual] <- y_d_cont[T_d == 0]
  y_cont_mat[3,1:n3] <- y_d_cont_rwd
  
  return(y_cont_mat)
}


for(iter in seq(0,99)){
#iter <- 0
  
  set.seed(12345+iter)
  
  n_RCT <- 300
  r <- 2
  n1 <- (r/(r+1))*n_RCT
  n2 <- (1/(r+1))*n_RCT
  n3 <- n_RCT
  
  p <- 3
  norm_weight <- c(0.33, 0.34, 0.33)
  norm_weight_RWD <- norm_weight
  
  # True clusters
  c1_mean <- rep(-2, p)
  c2_mean <- rep(2, p)
  allc_sig2 <- 1 
  
  cov_d <- matrix(rep(NA, (n1+n2)*p), ncol = p)
  c_d <- rep(NA, (n1+n2))
  T_d <- rep(NA, (n1+n2))
  
  for(i in 1:(n1+n2)){
    # Generate cluster membership
    c_d[i] <- sample(c(1,2,3), 1, replace = T, prob = norm_weight)
    # Generate treatment arm
    T_d[i] <- rbern(1, (r/(r+1)))
    # Generate covariates
    if(c_d[i] == 1){
      # In first cluster
      x2 <- rnorm(1,mean = -3,sd = 1)
      cov_d[i,] <- c(rnorm(1,mean = 2,sd = 1), x2, x2*m+rnorm(1))
    }else if(c_d[i] == 2){
      # In second cluster
      x2 <- rnorm(1,mean = 0,sd = 1)
      cov_d[i,] <- c(rnorm(1,mean = 0,sd = 1), x2, x2*m+rnorm(1))
    }else{
      # In third cluster
      x2 <- rnorm(1,mean = 3,sd = 1)
      cov_d[i,] <- c(rnorm(1,mean = -2,sd = 1), x2, x2*m+rnorm(1))
    }
  }
  n1_actual <- sum(T_d == 1)
  n2_actual <- sum(T_d == 0)
  
  cov_d_rwd <- matrix(rep(NA, n3*p), ncol = p)
  c_d_rwd <- rep(NA, n3)
  for(i in 1:n3){
    # Generate cluster membership
    c_d_rwd[i] <- sample(c(1,2,3), 1, replace = T, prob = norm_weight_RWD)
    # Generate covariates
    if(c_d_rwd[i] == 1){
      # In first cluster, common with G1/G2
      x2 <- rnorm(1,mean = -3,sd = 1)
      cov_d_rwd[i,] <- c(rnorm(1,mean = 2,sd = 1), x2, x2*m+rnorm(1))
    }else if(c_d_rwd[i] == 2){
      # In second cluster, common with G1/G2
      x2 <- rnorm(1,mean = 0,sd = 1)
      cov_d_rwd[i,] <- c(rnorm(1,mean = 0,sd = 1), x2, x2*m+rnorm(1))
    }else{
      # In third cluster, common with G1/G2
      x2 <- rnorm(1,mean = 3,sd = 1)
      cov_d_rwd[i,] <- c(rnorm(1,mean = -2,sd = 1), x2, x2*m+rnorm(1))
    }
  }
  
  # Input covariate and ground-truth clustering for PAM
  d_mat <- array(dim = c(3, max(c(n1_actual, n2_actual, n3)), p)) 
  c_mat <- matrix(rep(NA, max(c(n1_actual, n2_actual, n3))*3), nrow = 3) 
  
  d_mat[1,1:n1_actual,] <- cov_d[T_d == 1,]
  c_mat[1,1:n1_actual] <- c_d[T_d == 1]
  d_mat[2,1:n2_actual,] <- cov_d[T_d == 0,]
  c_mat[2,1:n2_actual] <- c_d[T_d == 0]
  d_mat[3,1:n3,] <- cov_d_rwd
  c_mat[3,1:n3] <- c_d_rwd
  
  saveRDS(list(d_mat = d_mat, c_mat = c_mat), file = paste0("./Sim_data/m100/Covariates_sim_",iter+1,".RData"))
  
  y_out <- gen_out(n3, n1_actual, n2_actual, 3, cov_d, T_d)
  saveRDS(y_out, file = paste0("./Sim_data/m100/Outcome_sim_",iter+1,".RData"))
  
}

