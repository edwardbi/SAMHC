rm(list = ls())

folder_dir <- "/Users/dehuabi/Desktop/" # Change this to where the folder is
work_dir <- paste0(folder_dir,"PAM_Design_Method")

setwd(work_dir)

# Load sources and libraries
library(binhf)
library(mclust)
library(mcclust)
library(mcclust.ext)
library(Rlab)
library(MASS)
library(boot)
library(psrwe)
library(invgamma)

set.seed(12345)

save_dir = "./results/"

##################################################################
# Read in simulated data
##################################################################

# Read in the data
file_name <- "sample.csv" # Change this to input file name
input_dir <- paste0("./input/",file_name)
samp_data <- read.csv(input_dir)

samp_data$EASI16 <- as.numeric(samp_data$EASI16)
samp_data2 <- data.frame(samp_data$STUDYID, 
                         samp_data$AGE, samp_data$EASIbl, samp_data$EASI16,
                         samp_data$height, samp_data$weight, samp_data$BSAbl)

samp_data2 <- samp_data2[complete.cases(samp_data2),]
samp_data2$bmi <- samp_data2$samp_data.weight/((samp_data2$samp_data.height/100)^2)

n1 <- sum(samp_data2$samp_data.STUDYID == "M16-045")
n2 <- sum(samp_data2$samp_data.STUDYID == "M18-891")
n3 <- sum(samp_data2$samp_data.STUDYID == "M16-047")

# Randomly sample half of group 2 to remove
n2_deduct <- n2 - floor(n1/2)
rm_g2_idx <- n1 + sample(seq(1,n2), n2_deduct, replace = F)
samp_data2_rmg2idx <- samp_data2[-rm_g2_idx, ]
n2_new <- (n2-length(rm_g2_idx))

d_mat <- array(dim = c(3, max(c(n1, n2_new, n3)), 4))
d_mat[1,1:n1,] <- as.matrix(samp_data2_rmg2idx[samp_data2_rmg2idx$samp_data.STUDYID == "M16-045",c(2:3,7:8)])
d_mat[2,1:n2_new,] <- as.matrix(samp_data2_rmg2idx[samp_data2_rmg2idx$samp_data.STUDYID == "M18-891",c(2:3,7:8)])
d_mat[3,1:n3,] <- as.matrix(samp_data2_rmg2idx[samp_data2_rmg2idx$samp_data.STUDYID == "M16-047",c(2:3,7:8)])

# Null outcome
easi75_n1 <- (samp_data2_rmg2idx$samp_data.EASIbl[1:n1] - samp_data2_rmg2idx$samp_data.EASI16[1:n1])/
  samp_data2_rmg2idx$samp_data.EASIbl[1:n1]
n1_y_null <- rep(0,n1)
n1_y_null[which(easi75_n1 >= 0.75)] <- 1
easi75_n2 <- (samp_data2_rmg2idx$samp_data.EASIbl[(n1+1):(n1+n2_new)] - samp_data2_rmg2idx$samp_data.EASI16[(n1+1):(n1+n2_new)])/
  samp_data2_rmg2idx$samp_data.EASIbl[(n1+1):(n1+n2_new)]
n2_y_null <- rep(0,n2_new)
n2_y_null[which(easi75_n2 >= 0.75)] <- 1
easi75_n3 <- (samp_data2_rmg2idx$samp_data.EASIbl[(n1+n2_new+1):(n1+n2_new+n3)] - samp_data2_rmg2idx$samp_data.EASI16[(n1+n2_new+1):(n1+n2_new+n3)])/
  samp_data2_rmg2idx$samp_data.EASIbl[(n1+n2_new+1):(n1+n2_new+n3)]
n3_y_null <- rep(0,n3)
n3_y_null[which(easi75_n3 >= 0.75)] <- 1
y_mat_null <- matrix(rep(NA, max(c(n1, n2_new, n3))*3), nrow = 3)
y_mat_null[1,1:n1] <- n1_y_null
y_mat_null[2,1:n2_new] <- n2_y_null
y_mat_null[3,1:n3] <- n3_y_null

# Generate alternative outcome
# Fit a regression with covariates and n1_y_null
fitted_model <- glm(n1_y_null ~ samp_data2_rmg2idx$samp_data.AGE[1:n1] + samp_data2_rmg2idx$samp_data.EASIbl[1:n1] +
                      samp_data2_rmg2idx$samp_data.BSAbl[1:n1] + samp_data2_rmg2idx$bmi[1:n1], family = "binomial")
model_y1_coeff <- as.vector(fitted_model$coefficients)
thres_val <- 0.8
delta <- NA
delta_seq <- seq(0, 10, by = 0.01)
n1_y_alt <- rep(0,n1)
for(delta_extra_idx in 1:length(delta_seq)){
  model_y1_coeff[1] <- model_y1_coeff[1] + delta_seq[delta_extra_idx]
  n1_y_alt_temp <- rep(0,n1)
  for(i in 1:n1){
    y_cov <- d_mat[1,i,]
    logit_val <- sum(y_cov*(model_y1_coeff[2:length(model_y1_coeff)])) + model_y1_coeff[1]
    n1_y_alt_temp[i] <- rbern(1, inv.logit(logit_val))
  }
  if((sum(n1_y_alt_temp)/length(n1_y_alt_temp)) > thres_val){
    n1_y_alt <- n1_y_alt_temp
    delta <- delta_seq[delta_extra_idx]
    break
  }
}

y_mat_alt <- matrix(rep(NA, max(c(n1, n2_new, n3))*3), nrow = 3)
y_mat_alt[1,1:n1] <- n1_y_alt
y_mat_alt[2,1:n2_new] <- n2_y_null
y_mat_alt[3,1:n3] <- n3_y_null

cov1_raw <- as.vector(t(d_mat[,,1]))
cov1_raw <- cov1_raw[complete.cases(cov1_raw)]
cov2_raw <- as.vector(t(d_mat[,,2]))
cov2_raw <- cov2_raw[complete.cases(cov2_raw)]
cov3_raw <- as.vector(t(d_mat[,,3]))
cov3_raw <- cov3_raw[complete.cases(cov3_raw)]
cov4_raw <- as.vector(t(d_mat[,,4]))
cov4_raw <- cov4_raw[complete.cases(cov4_raw)]

outcome_null <- as.vector(t(y_mat_null))
outcome_null <- outcome_null[complete.cases(outcome_null)]

outcome_alt <- as.vector(t(y_mat_alt))
outcome_alt <- outcome_alt[complete.cases(outcome_alt)]

dat_df <- data.frame(
  group = c(rep(1,n1), rep(1,n2_new), rep(0,n3)),
  arm = c(rep(1,n1), rep(0,n2_new), rep(0,n3)),
  cov1 = cov1_raw,
  cov2 = cov2_raw,
  cov3 = cov3_raw,
  cov4 = cov4_raw,
  outcome_null = outcome_null,
  outcome_alt = outcome_alt
)

# Use default LR model
ps_fit_default <- psrwe_est(dat_df, v_covs = paste("cov", 1:4, sep = ""), 
                            v_grp = "group", cur_grp_level = 1,
                            v_arm = "arm", ctl_arm_level = 0, nstrata = 5)
# composite likelihood
ps_default_bor_rct <- psrwe_borrow(ps_fit_default, total_borrow = n1 - n2_new, method = "distance")
trteff_default_null <- psrwe_compl(ps_default_bor_rct, v_outcome = "outcome_null", outcome_type = "binary")
trteff_default_alt <- psrwe_compl(ps_default_bor_rct, v_outcome = "outcome_alt", outcome_type = "binary")

# Use default RF model
ps_fit_default_rf <- psrwe_est(dat_df, v_covs = paste("cov", 1:4, sep = ""), ps_method = "randomforest", 
                               v_grp = "group", cur_grp_level = 1,
                               v_arm = "arm", ctl_arm_level = 0, nstrata = 5)
# composite likelihood
ps_default_rf_bor_rct <- psrwe_borrow(ps_fit_default_rf, total_borrow = n1 - n2_new, method = "distance")
trteff_default_rf_null <- psrwe_compl(ps_default_rf_bor_rct, v_outcome = "outcome_null", outcome_type = "binary")
trteff_default_rf_alt <- psrwe_compl(ps_default_rf_bor_rct, v_outcome = "outcome_alt", outcome_type = "binary")

# Model selection with GLM
full_glm_model <- glm(group ~ cov1 + cov2 + cov3 + cov4 + cov1^2 + cov2^2 + cov3^2 + cov4^2 +
                        cov1*cov2 + cov1*cov3 + cov1*cov4 + cov2*cov3 + cov2*cov4 + cov3*cov4 +
                        cov1*cov2*cov3 + cov1*cov2*cov4 + cov1*cov3*cov4 + cov2*cov3*cov4 + cov1*cov2*cov3*cov4, 
                      data = dat_df, family = "binomial")
model_sel_result <- stepAIC(full_glm_model, direction = 'backward', trace = 0)
model_sel <- model_sel_result$formula
ps_fit_model_sel <- psrwe_est(dat_df, ps_fml = as.formula(model_sel), v_covs = paste("cov", 1:4, sep = ""), 
                              v_grp = "group", cur_grp_level = 1,
                              v_arm = "arm", ctl_arm_level = 0, nstrata = 5)
ps_ms_bor_rct <- psrwe_borrow(ps_fit_model_sel, total_borrow = n1 - n2_new, method = "distance")
trteff_ms_null <- psrwe_compl(ps_ms_bor_rct, v_outcome = "outcome_null", outcome_type = "continuous")
trteff_ms_alt <- psrwe_compl(ps_ms_bor_rct, v_outcome = "outcome_alt", outcome_type = "binary")

save_out_dir <- paste0(save_dir,"PSCL_method_results.RData")
saveRDS(list(ps_fit_default = ps_fit_default, ps_default_bor_rct = ps_default_bor_rct, 
             trteff_default_null = trteff_default_null, trteff_default_alt = trteff_default_alt,
             ps_fit_default_rf = ps_fit_default_rf, ps_default_rf_bor_rct = ps_default_rf_bor_rct,
             trteff_default_rf_null = trteff_default_rf_null, trteff_default_rf_alt = trteff_default_rf_alt,
             formula = model_sel, ps_fit_model_sel = ps_fit_model_sel, ps_ms_bor_rct = ps_ms_bor_rct, 
             trteff_ms_null = trteff_ms_null, trteff_ms_alt = trteff_ms_alt), file = save_out_dir)





