library(binhf)
source("./codes/MSAM1_slice.R")
source("./codes/helper_funcs.R")
library(mclust)
library(mcclust)
library(mcclust.ext)
library(Rlab)
library(MASS)
library(boot)
library(invgamma)
library(psrwe)

get_input_data <- function(cov_data_dir, out_data_dir, sc, sim_n, delta){
  cov_input <- paste0(cov_data_dir,"Covariates_sim_sc",sc,"_",sim_n,".RData")
  cov_data <- readRDS(cov_input)
  out_input <- paste0(out_data_dir,"Outcome_sim_sc",sc,"_",sim_n,"_",delta,".RData")
  out_data <- readRDS(out_input)
  
  d_mat <- cov_data$d_mat
  c_mat <- cov_data$c_mat
  y_cont_mat <- out_data$y_cont
  y_bin_mat <- out_data$y_bin
  
  p <- length(d_mat[1,1,])
  r <- 2
  
  n1_actual <- length(c_mat[1,complete.cases(c_mat[1,])])
  n2_actual <- length(c_mat[2,complete.cases(c_mat[2,])])
  n3 <- length(c_mat[3,complete.cases(c_mat[3,])])
  
  cov1_raw <- as.vector(t(d_mat[,,1]))
  cov1_raw <- cov1_raw[complete.cases(cov1_raw)]
  cov2_raw <- as.vector(t(d_mat[,,2]))
  cov2_raw <- cov2_raw[complete.cases(cov2_raw)]
  cov3_raw <- as.vector(t(d_mat[,,3]))
  cov3_raw <- cov3_raw[complete.cases(cov3_raw)]
  
  true_c_vec <- as.vector(t(c_mat))
  true_c_vec <- true_c_vec[complete.cases(true_c_vec)]
  
  outcome_cont <- as.vector(t(y_cont_mat))
  outcome_cont <- outcome_cont[complete.cases(outcome_cont)]
  
  outcome_bi <- as.vector(t(y_bin_mat))
  outcome_bi <- outcome_bi[complete.cases(outcome_bi)]
  
  dat_df <- data.frame(
    group = c(rep(1,n1_actual), rep(1,n2_actual), rep(0,n3)),
    arm = c(rep(1,n1_actual), rep(0,n2_actual), rep(0,n3)),
    cov1 = cov1_raw,
    cov2 = cov2_raw,
    cov3 = cov3_raw,
    true_c = true_c_vec,
    outcome_cont = outcome_cont,
    outcome_bi = outcome_bi
  )
  return(list(df = dat_df, n1 = n1_actual, n2 = n2_actual))
}

compute_PSCL <- function(dat_df, n1_actual, n2_actual, delta){
  
  # Use default LR model
  ps_fit_default <- psrwe_est(dat_df, v_covs = paste("cov", 1:3, sep = ""), 
                              v_grp = "group", cur_grp_level = 1,
                              v_arm = "arm", ctl_arm_level = 0, nstrata = 5)
  ps_default_bor_rct <- psrwe_borrow(ps_fit_default, total_borrow = n1_actual - n2_actual, method = "distance")
  trteff_default_cont <- psrwe_compl(ps_default_bor_rct, v_outcome = "outcome_cont", outcome_type = "continuous")
  trteff_default_bin <- psrwe_compl(ps_default_bor_rct, v_outcome = "outcome_bi", outcome_type = "binary")
  #saveRDS(list(ps = ps_fit_default, borrow = ps_default_bor_rct, 
  #             trt_cont = trteff_default_cont, trt_bin = trteff_default_bin), file = save_PSCL_LR_TE_dir)
  
  # Use default RF model
  ps_fit_default_rf <- psrwe_est(dat_df, v_covs = paste("cov", 1:3, sep = ""), 
                                 ps_method = "randomforest", v_grp = "group", cur_grp_level = 1,
                                 v_arm = "arm", ctl_arm_level = 0, nstrata = 5)
  ps_default_rf_bor_rct <- psrwe_borrow(ps_fit_default_rf, total_borrow = n1_actual - n2_actual, method = "distance")
  trteff_default_rf_cont <- psrwe_compl(ps_default_rf_bor_rct, v_outcome = "outcome_cont", outcome_type = "continuous")
  trteff_default_rf_bin <- psrwe_compl(ps_default_rf_bor_rct, v_outcome = "outcome_bi", outcome_type = "binary")
  #saveRDS(list(ps = ps_fit_default_rf, borrow = ps_default_rf_bor_rct, 
  #             trt_cont = trteff_default_rf_cont, trt_bin = trteff_default_rf_bin), file = save_PSCL_RF_TE_dir)
  
  # Use model selection LR model
  full_glm_model <- glm(group ~ cov1 + cov2 + cov3 + cov1^2 + cov2^2 + cov3^2 +
                          cov1*cov2 + cov1*cov3 + cov2*cov3, data = dat_df, family = "binomial")
  model_sel_result <- stepAIC(full_glm_model, direction = 'backward', trace = 0)
  model_sel <- model_sel_result$formula
  ps_fit_model_sel <- psrwe_est(dat_df, ps_fml = as.formula(model_sel),
                                v_covs = paste("cov", 1:3, sep = ""), 
                                v_grp = "group", cur_grp_level = 1,
                                v_arm = "arm", ctl_arm_level = 0, nstrata = 5)
  ps_ms_bor_rct <- psrwe_borrow(ps_fit_model_sel, total_borrow = n1_actual - n2_actual, method = "distance")
  trteff_ms_cont <- psrwe_compl(ps_ms_bor_rct, v_outcome = "outcome_cont", outcome_type = "continuous")
  trteff_ms_bin <- psrwe_compl(ps_ms_bor_rct, v_outcome = "outcome_bi", outcome_type = "binary")
  #saveRDS(list(formula = model_sel, ps = ps_fit_model_sel, borrow = ps_ms_bor_rct, 
  #             trt_cont = trteff_ms_cont, trt_bin = trteff_ms_bin), file = save_PSCL_MS_TE_dir)
  
  trt_eff_cont <- c(trteff_default_cont$Effect$Overall_Estimate$Mean, 
                    trteff_default_rf_cont$Effect$Overall_Estimate$Mean,
                    trteff_ms_cont$Effect$Overall_Estimate$Mean)
  trt_eff_cont_sd <- c(trteff_default_cont$Effect$Overall_Estimate$StdErr, 
                       trteff_default_rf_cont$Effect$Overall_Estimate$StdErr,
                       trteff_ms_cont$Effect$Overall_Estimate$StdErr)
  bias_cont <- c(trt_eff_cont - rep(delta, 3))
  mse_cont <- bias_cont^2
  
  cont_eff <- list(trt = trt_eff_cont, sd = trt_eff_cont_sd, bias = bias_cont, mse = mse_cont)
  
  trt_eff_bin <- c(trteff_default_bin$Effect$Overall_Estimate$Mean, 
                   trteff_default_rf_bin$Effect$Overall_Estimate$Mean,
                   trteff_ms_bin$Effect$Overall_Estimate$Mean)
  trt_eff_bin_sd <- c(trteff_default_bin$Effect$Overall_Estimate$StdErr, 
                      trteff_default_rf_bin$Effect$Overall_Estimate$StdErr,
                      trteff_ms_bin$Effect$Overall_Estimate$StdErr)
  bias_bin <- c(trt_eff_bin - c(0,9.42,15.67,19.52)/100)
  mse_bin <- bias_bin^2
  
  bin_eff <- list(trt = trt_eff_bin, sd = trt_eff_bin_sd, bias = bias_bin, mse = mse_bin)
  
  return(list(cont = cont_eff, bin = bin_eff))
}

sc <- 3

## --- build file list ---
sim_id <- seq(1,1000)

cov_data_dir <- paste0("./data/Sc",sc,"_300/")
out_data_dir <- paste0("./data/Sc",sc,"_300/")
save_dir <- paste0("./results/Sc",sc,"_300/")

for(sim_n in sim_id){
  print(sim_n)
  dat_df_0 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 0)
  PSCL_0 <- compute_PSCL(dat_df_0$df, dat_df_0$n1, dat_df_0$n2, 0)
  
  dat_df_1 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 1)
  PSCL_1 <- compute_PSCL(dat_df_1$df, dat_df_1$n1, dat_df_1$n2, 1)
  
  dat_df_2 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 2)
  PSCL_2 <- compute_PSCL(dat_df_2$df, dat_df_2$n1, dat_df_2$n2, 2)
  
  dat_df_3 <- get_input_data(cov_data_dir, out_data_dir, sc, sim_n, 3)
  PSCL_3 <- compute_PSCL(dat_df_3$df, dat_df_3$n1, dat_df_3$n2, 3)
  
  PSCL_result <- list(PSCL0 = PSCL_0, PSCL1 = PSCL_1, PSCL2 = PSCL_2, PSCL3 = PSCL_3)
  save_PSCL_dir <- paste0(save_dir,"PSCL_sim_sc",sc,"_",sim_n,".RData")
  saveRDS(PSCL_result, file = save_PSCL_dir)
  
}


