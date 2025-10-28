rm(list = ls())

##################################################################
# Command line arguments
##################################################################

args = commandArgs(TRUE)

if(length(args) != 3){
  stop("Usage: SAM1_m1_script.R <seed_n> <work_dir>")
}

work_dir = args[1]
sc = args[2]
cov_data_dir = paste0("./data/Sc",sc,"_300/")
out_data_dir = paste0("./data/Sc",sc,"_300/")
save_dir = paste0("./results/Sc",sc,"_300/")
sim_n = args[3]

setwd(work_dir)
set.seed(1001)

##################################################################
# Pcakages and Sources
##################################################################


library(binhf)
library(mclust)
library(mcclust)
library(mcclust.ext)
library(Rlab)
library(MASS)
library(boot)
library(invgamma)

source("./codes/MSAM1_slice.R")
source("./codes/helper_funcs.R")

##################################################################
# Read in simulated data
##################################################################

cov_input <- paste0(cov_data_dir,"Covariates_sim_sc",sc,"_",sim_n,".RData")
cov_data <- readRDS(cov_input)

out_input_0 <- paste0(out_data_dir,"Outcome_sim_sc",sc,"_",sim_n,"_0.RData")
out_input_1 <- paste0(out_data_dir,"Outcome_sim_sc",sc,"_",sim_n,"_1.RData")
out_input_2 <- paste0(out_data_dir,"Outcome_sim_sc",sc,"_",sim_n,"_2.RData")
out_input_3 <- paste0(out_data_dir,"Outcome_sim_sc",sc,"_",sim_n,"_3.RData")

out_data_0 <- readRDS(out_input_0)
out_data_1 <- readRDS(out_input_1)
out_data_2 <- readRDS(out_input_2)
out_data_3 <- readRDS(out_input_3)

d_mat <- cov_data$d_mat
c_mat <- cov_data$c_mat
y_cont_mat_0 <- out_data_0$y_cont
y_bin_mat_0 <- out_data_0$y_bin
y_cont_mat_1 <- out_data_1$y_cont
y_bin_mat_1 <- out_data_1$y_bin
y_cont_mat_2 <- out_data_2$y_cont
y_bin_mat_2 <- out_data_2$y_bin
y_cont_mat_3 <- out_data_3$y_cont
y_bin_mat_3 <- out_data_3$y_bin

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

outcome_cont_0 <- as.vector(t(y_cont_mat_0))
outcome_cont_0 <- outcome_cont_0[complete.cases(outcome_cont_0)]
outcome_bi_0 <- as.vector(t(y_bin_mat_0))
outcome_bi_0 <- outcome_bi_0[complete.cases(outcome_bi_0)]

outcome_cont_1 <- as.vector(t(y_cont_mat_1))
outcome_cont_1 <- outcome_cont_1[complete.cases(outcome_cont_1)]
outcome_bi_1 <- as.vector(t(y_bin_mat_1))
outcome_bi_1 <- outcome_bi_1[complete.cases(outcome_bi_1)]

outcome_cont_2 <- as.vector(t(y_cont_mat_2))
outcome_cont_2 <- outcome_cont_2[complete.cases(outcome_cont_2)]
outcome_bi_2 <- as.vector(t(y_bin_mat_2))
outcome_bi_2 <- outcome_bi_2[complete.cases(outcome_bi_2)]

outcome_cont_3 <- as.vector(t(y_cont_mat_3))
outcome_cont_3 <- outcome_cont_3[complete.cases(outcome_cont_3)]
outcome_bi_3 <- as.vector(t(y_bin_mat_3))
outcome_bi_3 <- outcome_bi_3[complete.cases(outcome_bi_3)]

# print(c(length(outcome_bi_0), length(outcome_bi_1), length(outcome_bi_2), length(outcome_bi_3)))

dat_df <- data.frame(
  group = c(rep(1,n1_actual), rep(1,n2_actual), rep(0,n3)),
  arm = c(rep(1,n1_actual), rep(0,n2_actual), rep(0,n3)),
  cov1 = cov1_raw,
  cov2 = cov2_raw,
  cov3 = cov3_raw,
  true_c = true_c_vec,
  outcome_cont_0 = outcome_cont_0,
  outcome_bi_0 = outcome_bi_0,
  outcome_cont_1 = outcome_cont_1,
  outcome_bi_1 = outcome_bi_1,
  outcome_cont_2 = outcome_cont_2,
  outcome_bi_2 = outcome_bi_2,
  outcome_cont_3 = outcome_cont_3,
  outcome_bi_3 = outcome_bi_3
)

se_baseline_bi <- function(dat_df, n1_actual, n2_actual){
  
  a <- sum(dat_df$outcome_bi[dat_df$arm == 1])
  #b <- n1_actual - a
  #apb <- n1_actual
  c <- sum(dat_df$outcome_bi[dat_df$group == 1 & dat_df$arm == 0])
  #d <- n2_actual - c
  #cpd <- n2_actual
  #apc <- a+c
  #bpd <- b+d
  #N <- n1_actual+n2_actual
  #se <- sqrt(apc*bpd/((N^2)*apb))
  
  p1 <- a/n1_actual
  p2 <- c/n2_actual
  
  se <- sqrt(p1*(1-p1)/n1_actual + p2*(1-p2)/n2_actual)
  #print(se)
  wald_ts <- (p1 - p2)/se
  p_val <- pnorm(wald_ts, lower.tail = F)
  
  return(list(se = se, p_val = p_val))
  
}

comp_cov_prob <- function(trt_eff_MCMC, MCMC_itr, true_v, creb_int = 0.95){
  lower_cut <- (1 - creb_int)/2
  upper_cut <- (1 + creb_int)/2
  sorted_MCMC <- sort(trt_eff_MCMC, decreasing = F)
  sorted_MCMC <- na.omit(sorted_MCMC)
  lower_cred_bd <- sorted_MCMC[floor(length(sorted_MCMC)*lower_cut)]
  upper_cred_bd <- sorted_MCMC[ceiling(length(sorted_MCMC)*upper_cut)]
  if(true_v >= lower_cred_bd && true_v <= upper_cred_bd){
    return(list(cp = 1, lower_cbd = lower_cred_bd, upper_cbd = upper_cred_bd))
  }else{
    return(list(cp = 0, lower_cbd = lower_cred_bd, upper_cbd = upper_cred_bd))
  }
}

comp_cov_prob_om <- function(mean, sd, true_v, cv = 1.96){
  lower_ci <- mean - sd*cv
  upper_ci <- mean + sd*cv
  if(true_v >= lower_ci && true_v <= upper_ci){
    return(list(cp = 1, lower_cbd = lower_ci, upper_cbd = upper_ci))
  }else{
    return(list(cp = 0, lower_cbd = lower_ci, upper_cbd = upper_ci))
  }
}

power_comp <- function(trt_eff_MCMC, level = 0.95){
  
  sorted_MCMC <- sort(trt_eff_MCMC, decreasing = F)
  sorted_MCMC <- na.omit(sorted_MCMC)
  prob_g0 <- length(which(sorted_MCMC > 0))/length(sorted_MCMC)
  if(prob_g0 >= level){
    return(list(dec = 1, prob = prob_g0))
  }else{
    return(list(dec = 0, prob = prob_g0))
  }
  
}

power_comp_om <- function(mean, sd, level = 0.05){
  
  z_score <- mean/sd
  p_val <- pnorm(z_score, lower.tail = F)
  if(p_val <= level){
    return(list(dec = 1, p_val = p_val))
  }else{
    return(list(dec = 0, p_val = p_val))
  }
  
}

truth_eff_cont <- c(0, 1, 2, 3)
truth_eff_bin <- rep(NA, 4)
for(i in 1:4){
  if(i == 1){
    truth_eff_bin[i] <- 0
  }else{
    truth_eff_bin[i] <- (0.3*inv.logit(6+truth_eff_cont[i]) + 0.4*inv.logit(truth_eff_cont[i]) + 0.3*inv.logit(-6+truth_eff_cont[i])) -
      (0.3*inv.logit(6) + 0.4*inv.logit(0) + 0.3*inv.logit(-6))
  }
}

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
#post_result_file <- paste0("results/Sc",sc,"_300/PAM_design_MCMC_sim_sc",sc,"_",sim_n,".RData")
#post_result <- readRDS(post_result_file)
post_result <- MCMC_SAM1(d_mat, mu0, sig0, a, b, kappa0, nu0, a1, b1, a2, b2, 
                         0.5, 0.3, 0.5, 10, M_iter, burnin)
#save_PAM_MCMC_dir <- paste0(save_dir,"PAM_design_MCMC_sim_sc",sc,"_",sim_n,".RData")
#saveRDS(post_result, file = save_PAM_MCMC_dir)

# Hyperparameters for PAM design
A <- n1_actual - n2_actual
a0 <- 0.5
b0 <- 0.5
a0_y <- 3
b0_y <- 3
mu0_y <- 0
nu0_y <- 0.1

PAM_result_0 <- PAM_design_pp(d_mat, y_cont_mat_0, y_bin_mat_0, post_result, n1_actual, n2_actual, n3, 
                              r, M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)
#save_PAM_TE_dir_0 <- paste0(save_dir,"PAM_design_TE_sim_sc",sc,"_",sim_n,"_0.RData")
#saveRDS(PAM_result_0, file = save_PAM_TE_dir_0)

PAM_result_1 <- PAM_design_pp(d_mat, y_cont_mat_1, y_bin_mat_1, post_result, n1_actual, n2_actual, n3, 
                              r, M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)
#save_PAM_TE_dir_1 <- paste0(save_dir,"PAM_design_TE_sim_sc",sc,"_",sim_n,"_1.RData")
#saveRDS(PAM_result_1, file = save_PAM_TE_dir_1)

PAM_result_2 <- PAM_design_pp(d_mat, y_cont_mat_2, y_bin_mat_2, post_result, n1_actual, n2_actual, n3, 
                              r, M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)
#save_PAM_TE_dir_2 <- paste0(save_dir,"PAM_design_TE_sim_sc",sc,"_",sim_n,"_2.RData")
#saveRDS(PAM_result_2, file = save_PAM_TE_dir_2)

PAM_result_3 <- PAM_design_pp(d_mat, y_cont_mat_3, y_bin_mat_3, post_result, n1_actual, n2_actual, n3, 
                              r, M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)
#save_PAM_TE_dir_3 <- paste0(save_dir,"PAM_design_TE_sim_sc",sc,"_",sim_n,"_3.RData")
#saveRDS(PAM_result_3, file = save_PAM_TE_dir_3)

PAM_trt <- c(mean(PAM_result_0$trt_eff_cont[2000:5000], na.rm = T),
             mean(PAM_result_1$trt_eff_cont[2000:5000], na.rm = T),
             mean(PAM_result_2$trt_eff_cont[2000:5000], na.rm = T),
             mean(PAM_result_3$trt_eff_cont[2000:5000], na.rm = T))

PAM_trt_sd <- c(sd(PAM_result_0$trt_eff_cont[2000:5000], na.rm = T),
                sd(PAM_result_1$trt_eff_cont[2000:5000], na.rm = T),
                sd(PAM_result_2$trt_eff_cont[2000:5000], na.rm = T),
                sd(PAM_result_3$trt_eff_cont[2000:5000], na.rm = T))

cov_prob <- c(comp_cov_prob(PAM_result_0$trt_eff_cont[2001:5000], 3000, 0)$cp, 
              comp_cov_prob(PAM_result_1$trt_eff_cont[2001:5000], 3000, 1)$cp,
              comp_cov_prob(PAM_result_2$trt_eff_cont[2001:5000], 3000, 2)$cp, 
              comp_cov_prob(PAM_result_3$trt_eff_cont[2001:5000], 3000, 3)$cp)

power_val <- c(power_comp(PAM_result_0$trt_eff_cont[2001:5000])$dec, 
               power_comp(PAM_result_1$trt_eff_cont[2001:5000])$dec,
               power_comp(PAM_result_2$trt_eff_cont[2001:5000])$dec, 
               power_comp(PAM_result_3$trt_eff_cont[2001:5000])$dec)

save_PAM_data_cont <- list(trt_eff = PAM_trt, sds_val = PAM_trt_sd,
                           bias = PAM_trt - truth_eff_cont, 
                           mse = (PAM_trt - truth_eff_cont)^2,
                           cov_prob = cov_prob, power = power_val)

PAM_trt_bin <- c(mean(PAM_result_0$trt_eff_bin[2000:5000], na.rm = T),
                 mean(PAM_result_1$trt_eff_bin[2000:5000], na.rm = T),
                 mean(PAM_result_2$trt_eff_bin[2000:5000], na.rm = T),
                 mean(PAM_result_3$trt_eff_bin[2000:5000], na.rm = T))

PAM_trt_sd_bin <- c(sd(PAM_result_0$trt_eff_bin[2000:5000], na.rm = T),
                    sd(PAM_result_1$trt_eff_bin[2000:5000], na.rm = T),
                    sd(PAM_result_2$trt_eff_bin[2000:5000], na.rm = T),
                    sd(PAM_result_3$trt_eff_bin[2000:5000], na.rm = T))

cov_prob_bin <- c(comp_cov_prob(PAM_result_0$trt_eff_bin[2001:5000], 3000, truth_eff_bin[1])$cp, 
                  comp_cov_prob(PAM_result_1$trt_eff_bin[2001:5000], 3000, truth_eff_bin[2])$cp,
                  comp_cov_prob(PAM_result_2$trt_eff_bin[2001:5000], 3000, truth_eff_bin[3])$cp, 
                  comp_cov_prob(PAM_result_3$trt_eff_bin[2001:5000], 3000, truth_eff_bin[4])$cp)

power_val_bin <- c(power_comp(PAM_result_0$trt_eff_bin[2001:5000])$dec, 
                   power_comp(PAM_result_1$trt_eff_bin[2001:5000])$dec,
                   power_comp(PAM_result_2$trt_eff_bin[2001:5000])$dec, 
                   power_comp(PAM_result_3$trt_eff_bin[2001:5000])$dec)

save_PAM_data_bin <- list(trt_eff = PAM_trt_bin, sds_val = PAM_trt_sd_bin,
                          bias = PAM_trt_bin - truth_eff_bin, 
                          mse = (PAM_trt_bin - truth_eff_bin)^2,
                          cov_prob = cov_prob_bin, power = power_val_bin)

save_PAM_data <- list(cont = save_PAM_data_cont, 
                      bin = save_PAM_data_bin)

save_PAM_TE_dir <- paste0(save_dir,"PAM_design_TE_sim_sc",sc,"_",sim_n,".RData")
saveRDS(save_PAM_data, file = save_PAM_TE_dir)

##################################################################
# Baseline effect
##################################################################

bi_te_baseline <- function(y_bin_m){
  
  baseline_te_bin_trt <- glm(y_bin_m[1,] ~ 1, family = "binomial")
  te_bin_summary_trt <- summary(baseline_te_bin_trt)
  baseline_te_bin_ctl <- glm(y_bin_m[2,] ~ 1, family = "binomial")
  te_bin_summary_ctl <- summary(baseline_te_bin_ctl)
  te_bin_summary <- inv.logit(te_bin_summary_trt$coefficients[1,1]) - inv.logit(te_bin_summary_ctl$coefficients[1,1])
  p1 <- inv.logit(te_bin_summary_trt$coefficients[1,1])
  p2 <- inv.logit(te_bin_summary_ctl$coefficients[1,1])
  se <- sqrt(p1*(1-p1)/length(y_bin_m[1,complete.cases(y_bin_m[1,])]) + 
               p2*(1-p2)/length(y_bin_m[2,complete.cases(y_bin_m[2,])]))
  
  return(list(te_bin = te_bin_summary, se_bin = se))
}

save_baseline_TE_dir <- paste0(save_dir,"Baseline_TE_sim_sc",sc,"_",sim_n,".RData")

# Estimate baseline overall treatment effect
baseline_te_cont_0 <- lm(dat_df$outcome_cont_0[dat_df$group == 1] ~ dat_df$arm[dat_df$group == 1])
te_cont_summary_0 <- summary(baseline_te_cont_0)
baseline_te_cont_1 <- lm(dat_df$outcome_cont_1[dat_df$group == 1] ~ dat_df$arm[dat_df$group == 1])
te_cont_summary_1 <- summary(baseline_te_cont_1)
baseline_te_cont_2 <- lm(dat_df$outcome_cont_2[dat_df$group == 1] ~ dat_df$arm[dat_df$group == 1])
te_cont_summary_2 <- summary(baseline_te_cont_2)
baseline_te_cont_3 <- lm(dat_df$outcome_cont_3[dat_df$group == 1] ~ dat_df$arm[dat_df$group == 1])
te_cont_summary_3 <- summary(baseline_te_cont_3)

baseline_te_bin_0 <- bi_te_baseline(y_bin_mat_0)
baseline_te_bin_1 <- bi_te_baseline(y_bin_mat_1)
baseline_te_bin_2 <- bi_te_baseline(y_bin_mat_2)
baseline_te_bin_3 <- bi_te_baseline(y_bin_mat_3)

saveRDS(list(trt_cont_0 = te_cont_summary_0, trt_cont_1 = te_cont_summary_1, trt_cont_2 = te_cont_summary_2, trt_cont_3 = te_cont_summary_3,
             trt_bin_0 = baseline_te_bin_0, trt_bin_1 = baseline_te_bin_1, trt_bin_2 = baseline_te_bin_2, trt_bin_3 = baseline_te_bin_3), 
        file = save_baseline_TE_dir)

