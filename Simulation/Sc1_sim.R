rm(list = ls())

setwd("/Users/dehuabi/Desktop/PAM_Design")

library(Rlab)
library(MASS)
library(boot)

for(iter in seq(1,1000)){
  #iter <- 0
  set.seed(12345+iter)
  
  # ##################################################################
  # #################### Generate Covariates Data ####################
  # ##################################################################
  
  # Simulate three groups:
  # n1 = 200, n2 = 100, n3 = 300
  # G1, G2 same distribution (mixture of Normal)
  # G3 contains all clusters in G1/G2 with extra clusters
  
  # Generate current RCT
  # n_RCT refers to number of paitents in current RCT
  # n1 refers to treatment group
  # n2 refers to control group
  # r refers to randomization ratio: 
  #     Example: 2:1 randomization,
  #              r = 2
  # p refers to covariate dimension
  # norm_weight refers to ground-truth normal mixture weights
  # delta refers to the treatment effect
  n_RCT <- 450
  r <- 2
  n1 <- (r/(r+1))*n_RCT
  n2 <- (1/(r+1))*n_RCT
  p <- 3
  norm_weight <- c(0.3, 0.4, 0.3)
  n3 <- n_RCT
  norm_weight_RWD <- c(0.2, 0.3, 0.3, 0.2)
  
  # Truth clusters
  c1_mean <- rep(2, p) 
  c2_mean <- rep(0, p)
  c3_mean <- rep(-2,p) 
  c4_mean <- rep(-4,p) 
  allc_sig2 <- 1       
  
  # cov_d contains matrix of covariates in G1 and G2
  # cov_d_rwd contains matrix of covariates in G3
  # c_d contains vector of cluster membership in G1 and G2
  # c_d_rwd contains vector of cluster membership in G3
  # T_d contains randomization to treatment (grp 1) or ctl (grp 2)
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
      cov_d[i,] <- MASS::mvrnorm(n = 1, c1_mean, allc_sig2*diag(nrow = p))
    }else if(c_d[i] == 2){
      # In second cluster
      cov_d[i,] <- MASS::mvrnorm(n = 1, c2_mean, allc_sig2*diag(nrow = p))
    }else{
      # In third cluster
      cov_d[i,] <- MASS::mvrnorm(n = 1, c3_mean, allc_sig2*diag(nrow = p))
    }
  }
  n1_actual <- sum(T_d == 1)
  n2_actual <- sum(T_d == 0)
  
  cov_d_rwd <- matrix(rep(NA, n3*p), ncol = p)
  c_d_rwd <- rep(NA, n3)
  for(i in 1:n3){
    # Generate cluster membership
    c_d_rwd[i] <- sample(c(1,2,3,4), 1, replace = T, prob = norm_weight_RWD)
    # Generate covariates
    if(c_d_rwd[i] == 1){
      # In first cluster, common with G1/G2
      cov_d_rwd[i,] <- MASS::mvrnorm(n = 1, c1_mean, allc_sig2*diag(nrow = p))
    }else if(c_d_rwd[i] == 2){
      # In second cluster, common with G1/G2
      cov_d_rwd[i,] <- MASS::mvrnorm(n = 1, c2_mean, allc_sig2*diag(nrow = p))
    }else if(c_d_rwd[i] == 3){
      # In third cluster, common with G1/G2
      cov_d_rwd[i,] <- MASS::mvrnorm(n = 1, c3_mean, allc_sig2*diag(nrow = p))
    }else{
      # In fourth cluster, unique in G3
      cov_d_rwd[i,] <- MASS::mvrnorm(n = 1, c4_mean, allc_sig2*diag(nrow = p))
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
  
  #saveRDS(list(d_mat = d_mat, c_mat = c_mat), file = paste0("./Sim_data/Sc1_Data_e1/Covariates_sim_sc1_",iter+1,".RData"))
  saveRDS(list(d_mat = d_mat, c_mat = c_mat), file = paste0("./Sc1_",n_RCT,"/Covariates_sim_sc1_",iter,".RData"))
  
  
  # ##################################################################
  # ####################  Generate Outcomes Data  ####################
  # ##################################################################
  
  # Treatment effect
  gen_out <- function(n1, n2, n3, n1_actual, n2_actual, delta, cov_d, T_d){
    # y_d_cont, y_d_bin the continuous/binary outcome in G1 and G2
    # y_d_cont_rwd, y_d_bin_rwd the continuous/binary outcomes in G3
    y_d_cont <- rep(NA, (n1+n2))
    y_d_bin <- rep(NA, (n1+n2))
    for(i in 1:(n1+n2)){
      if(T_d[i] == 1){
        # Treatment arm
        y_d_cont[i] <- sum(cov_d[i,]) + delta + rnorm(1)
        y_bin_prob <- inv.logit(sum(cov_d[i,]) + delta)
        y_d_bin[i] <- rbern(1, y_bin_prob)
      }else{
        # Control arm
        y_d_cont[i] <- sum(cov_d[i,]) + rnorm(1)
        y_bin_prob <- inv.logit(sum(cov_d[i,]))
        y_d_bin[i] <- rbern(1, y_bin_prob)
      }
    }
    y_d_cont_rwd <- rep(NA, n3)
    y_d_bin_rwd <- rep(NA, n3)
    for(i in 1:n3){
      # Generate continuous and binary outcomes
      y_d_cont_rwd[i] <- sum(cov_d_rwd[i,]) + rnorm(1)
      y_bin_rwd_prob <- inv.logit(sum(cov_d_rwd[i,]))
      y_d_bin_rwd[i] <- rbern(1, y_bin_rwd_prob)
    }
    
    y_cont_mat <- matrix(rep(NA, max(c(n1_actual, n2_actual, n3))*3), nrow = 3)
    y_cont_mat[1,1:n1_actual] <- y_d_cont[T_d == 1]
    y_cont_mat[2,1:n2_actual] <- y_d_cont[T_d == 0]
    y_cont_mat[3,1:n3] <- y_d_cont_rwd
    
    y_bin_mat <- matrix(rep(NA, max(c(n1_actual, n2_actual, n3))*3), nrow = 3)
    y_bin_mat[1,1:n1_actual] <- y_d_bin[T_d == 1]
    y_bin_mat[2,1:n2_actual] <- y_d_bin[T_d == 0]
    y_bin_mat[3,1:n3] <- y_d_bin_rwd
    
    return(list(y_cont = y_cont_mat, y_bin = y_bin_mat))
    
  }
  
  y_out_0 <- gen_out(n1, n2, n3, n1_actual, n2_actual, 0, cov_d, T_d)
  saveRDS(y_out_0, file = paste0("./Sc1_",n_RCT,"/Outcome_sim_sc1_",iter,"_0.RData"))
  
  y_out_1 <- gen_out(n1, n2, n3, n1_actual, n2_actual, 1, cov_d, T_d)
  saveRDS(y_out_1, file = paste0("./Sc1_",n_RCT,"/Outcome_sim_sc1_",iter,"_1.RData"))
  
  y_out_2 <- gen_out(n1, n2, n3, n1_actual, n2_actual, 2, cov_d, T_d)
  saveRDS(y_out_2, file = paste0("./Sc1_",n_RCT,"/Outcome_sim_sc1_",iter,"_2.RData"))
  
  y_out_3 <- gen_out(n1, n2, n3, n1_actual, n2_actual, 3, cov_d, T_d)
  saveRDS(y_out_3, file = paste0("./Sc1_",n_RCT,"/Outcome_sim_sc1_",iter,"_3.RData"))
  
  
  
  #y_out <- gen_out(n1, n2, n3, n1_actual, n2_actual, 3, cov_d, T_d)
  #saveRDS(y_out, file = paste0("./Sim_data/Sc1_Data_e1/Outcome_sim_sc1_",iter+1,".RData"))
  
  #y_out_null <- gen_out(n1, n2, n3, n1_actual, n2_actual, 2, cov_d, T_d)
  #saveRDS(y_out_null, file = paste0("./Sc1_300/Outcome_sim_sc1_",iter,"_2.RData"))
  
  #y_out_null <- gen_out(n1, n2, n3, n1_actual, n2_actual, 0, cov_d, T_d)
  #saveRDS(y_out_null, file = paste0("./Sim_data/Sc1_Data_e1/Outcome_sim_sc1_",iter+1,"_null.RData"))
}