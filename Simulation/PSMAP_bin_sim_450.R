
source("./codes/PSMAP_Binary/PSMAP_Binary/ps-power_functions.R") 
source("./codes/PSMAP_Binary/PSMAP_Binary/ps-power_stan.R")
source("./codes/PSMAP_Binary/PSMAP_Binary/PSMAP_Binary_functions.R")

sc <- 3

get_input_data <- function(cov_data_dir, out_data_dir, sc, sim_n, delta){
  cov_input <- paste0(cov_data_dir,"Covariates_sim_sc",sc,"_",sim_n,".RData")
  cov_data <- readRDS(cov_input)
  out_input <- paste0(out_data_dir,"Outcome_sim_sc",sc,"_",sim_n,"_",delta,".RData")
  out_data <- readRDS(out_input)
  
  return(list(cov = cov_data, out = out_data))
}

PSMAP_result <- function(cov_data, outcome_data, delta){
  
  n_trt <- length((cov_data$c_mat[1,])[complete.cases(cov_data$c_mat[1,])])
  trt_data <- cov_data$d_mat[1,1:n_trt,]
  n_ctl <- length((cov_data$c_mat[2,])[complete.cases(cov_data$c_mat[2,])])
  ctl_data <- cov_data$d_mat[2,1:n_ctl,]
  his_data <- cov_data$d_mat[3,,]
  n_his_tot <- n_trt + n_ctl
  
  p <- length(cov_data$d_mat[1,1,])
  r <- 2
  
  out_data_trt <- outcome_data$y_bin[1,1:n_trt]
  ybar_trt <- sum(out_data_trt)
  out_data_ctl <- outcome_data$y_bin[2,1:n_ctl]
  ctl_y_vec <- matrix(rep(NA,n_ctl),ncol = 1)
  ctl_y_vec[,1] <- out_data_ctl
  ybar_ctl <- sum(out_data_ctl)
  out_data_his <- outcome_data$y_bin[3,]
  his_y_vec <- matrix(rep(NA,n_his_tot),ncol = 1)
  his_y_vec[,1] <- out_data_his
  ybar_his <- sum(out_data_his)
  
  # PS MAP
  regroup_data <- PS_regroup(his_data, his_y_vec, ctl_data, ctl_y_vec, p = p, S = 2)
  if(0 %in% regroup_data$rS){
    regroup_data$rS[which(regroup_data$rS == 0)] <- rep(1e-4,length(which(regroup_data$rS == 0))) # uncomment this for all sc3
    #print(i)                                                                                     # Use these two lines for other scs
    #next
  }
  
  PS_MAP_fit <- PS_MAP.fit2(ybar.trt = ybar_trt, n.trt = n_trt, tau.init = 1, target.ESS = n_trt - n_ctl, n.cur = regroup_data$n.cur, 
                            ybar.hist = regroup_data$Ybar.hist, n.hist = regroup_data$n.hist, overlap_coefficients = regroup_data$rS,
                            lim = c(0.001, 2), data.indiv = ctl_y_vec, alpha0 = 1, beta0 = 1, niter = 10000)
  
  PS_MAP_model <- PS_MAP_fit[[1]]
  PS_MAP_ESS <- PS_MAP_fit[[2]]
  
  results <- PS_MAP_model[[1]][,"effect"]
  
  trt_eff <- mean(results[1])
  mse_eff <- (mean(results[1]) - delta)^2
  
  return(list(trt = trt_eff, mse = mse_eff))
}

sim_id <- seq(1, 1000)

cov_data_dir <- paste0("./data/Sc",sc,"_450/")
out_data_dir <- paste0("./data/Sc",sc,"_450/")
save_dir <- paste0("./results/Sc",sc,"_450/")

true_tox <- c(0,9.42,15.67,19.52)/100

for(sim_n in sim_id){
  print(sim_n)
  
  dat_df_0 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 0)
  PSMAP_0 <- PSMAP_result(dat_df_0$cov, dat_df_0$out, true_tox[1])
  
  dat_df_1 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 1)
  PSMAP_1 <- PSMAP_result(dat_df_1$cov, dat_df_1$out, true_tox[2])
  
  dat_df_2 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 2)
  PSMAP_2 <- PSMAP_result(dat_df_2$cov, dat_df_2$out, true_tox[3])
  
  dat_df_3 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 3)
  PSMAP_3 <- PSMAP_result(dat_df_3$cov, dat_df_3$out, true_tox[4])
  
  PSMAP_bin <- list(PSM0 = PSMAP_0, PSM1 = PSMAP_1, 
                    PSM2 = PSMAP_2, PSM3 = PSMAP_3)
  
  save_PSMAP_bin_dir <- paste0(save_dir,"PSMAP_bin_sim_sc",sc,"_",sim_n,".RData")
  saveRDS(PSMAP_bin, file = save_PSMAP_bin_dir)
}



