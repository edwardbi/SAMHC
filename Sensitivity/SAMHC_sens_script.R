rm(list = ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

#if(length(args) != 3){
#  stop("Usage: SAM1_m1_script.R <seed_n> <work_dir>")
#}

work_dir = "/Users/dehuabi/Desktop/PAMHC_sensitivity/"#args[1]
folder = args[1]
cov_data_dir = paste0("./Sim_data/",folder,"/")
out_data_dir = paste0("./Sim_data/",folder,"/")
save_dir = paste0("./trt_results/",folder,"/")
idx_s = as.integer(args[2])
idx_e = as.integer(args[3])

setwd(work_dir)
set.seed(12345)

mcmc_dir <- paste0(work_dir,"Results/",folder)
mcmc_files <- list.files(path=mcmc_dir, pattern="sim_*")#, full.names=TRUE, recursive=FALSE)
mcmc_files <- mcmc_files[(length(mcmc_files)/2+1):length(mcmc_files)]
coupled_opt_files <- list.files(path=mcmc_dir, pattern="opt_cluster_sim_")#, full.names=TRUE, recursive=FALSE)


##################################################################
# Pcakages and Sources
##################################################################

library(binhf)
library(overlapping)
source("./MSAM1_slice.R")
source("./helper_funcs.R")
source("./SAMHC_sens_core.R")
library(mclust)
library(mcclust)
library(mcclust.ext)
library(Rlab)
library(MASS)
library(boot)
library(invgamma)
library(stringr)

##################################################################
# Read in simulated data
##################################################################

regexp <- "[[:digit:]]+"

for(idx in seq(idx_s,idx_e)){
  #idx <- 1
  sim_n <- as.integer(str_extract(mcmc_files[idx], regexp))
  cov_input <- paste0(cov_data_dir,"Covariates_sim_",sim_n,".RData")
  cov_data <- readRDS(cov_input)
  out_input <- paste0(out_data_dir,"Outcome_sim_",sim_n,".RData")
  out_data <- readRDS(out_input)
  
  d_mat <- cov_data$d_mat[,,c(1,2)]
  c_mat <- cov_data$c_mat
  y_cont_mat <- out_data
  
  p <- length(d_mat[1,1,])
  r <- 2
  
  n1_actual <- length(c_mat[1,complete.cases(c_mat[1,])])
  n2_actual <- length(c_mat[2,complete.cases(c_mat[2,])])
  n3 <- length(c_mat[3,complete.cases(c_mat[3,])])
  
  cov1_raw <- as.vector(t(d_mat[,,1]))
  cov1_raw <- cov1_raw[complete.cases(cov1_raw)]
  cov2_raw <- as.vector(t(d_mat[,,2]))
  cov2_raw <- cov2_raw[complete.cases(cov2_raw)]
  
  true_c_vec <- as.vector(t(c_mat))
  true_c_vec <- true_c_vec[complete.cases(true_c_vec)]
  
  outcome_cont <- as.vector(t(y_cont_mat))
  outcome_cont <- outcome_cont[complete.cases(outcome_cont)]
  
  dat_df <- data.frame(
    group = c(rep(1,n1_actual), rep(1,n2_actual), rep(0,n3)),
    arm = c(rep(1,n1_actual), rep(0,n2_actual), rep(0,n3)),
    cov1 = cov1_raw,
    cov2 = cov2_raw,
    true_c = true_c_vec,
    outcome_cont = outcome_cont
  )

  ##################################################################
  # PAM Result effect
  ##################################################################
  
  # Hyperparameters
  M_iter <- 5000
  burnin <- 5000
  
  mu0 <- rep(0, p)
  sig0 <- diag(nrow = p)
  
  a <- 0.5
  b <- 0.5
  
  kappa0 <- 0.1
  nu0 <- p
  
  a1 <- 3
  b1 <- 3
  
  a2 <- 3
  b2 <- 3
  
  # Step 1: Cluster treatment, control, and RWD with PAM
  post_result <- readRDS(paste0("./Results/",folder,"/",mcmc_files[idx]))
  clust_result <- readRDS(paste0("./Results/",folder,"/",coupled_opt_files[idx]))
  opt_clust <- clust_result$optclust$cl[1,]
  
  # Hyperparameters for PAM design
  A <- n1_actual - n2_actual
  a0 <- 0.5
  b0 <- 0.5
  a0_y <- 3
  b0_y <- 3
  mu0_y <- 0
  nu0_y <- 0.1
  
  PAM_result <- PAM_design_pp_y_sens(d_mat, y_cont_mat, post_result, n1_actual, n2_actual, n3, r,
                                     M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)
  save_PAM_TE_dir <- paste0(save_dir, "PAM_design_TE_sim_",folder,"_",sim_n,".RData")
  saveRDS(PAM_result, file = save_PAM_TE_dir)
  
  # Test-then-pool
  SAM_TTP_result <- SAM_TTP(d_mat, y_cont_mat, opt_clust, n1_actual, n2_actual, n3, r,
                            a0, b0, a0_y, b0_y, mu0_y, nu0_y)
  save_PAMTTE_TE_dir <- paste0(save_dir, "PAM_TTE_TE_sim_",folder,"_",sim_n,".RData")
  saveRDS(SAM_TTP_result, file = save_PAMTTE_TE_dir)
  
  # Baseline
  baseline_te_cont <- lm(dat_df$outcome_cont[dat_df$group == 1] ~ dat_df$arm[dat_df$group == 1])
  te_cont_summary <- summary(baseline_te_cont)
  trt_cont <- list(mean = te_cont_summary$coefficients[2,1], sd = te_cont_summary$coefficients[2,2])
  save_BL_TE_dir <- paste0(save_dir, "BL_TE_sim_",folder,"_",sim_n,".RData")
  saveRDS(trt_cont, file = save_BL_TE_dir)
  
  print(sim_n)
}

