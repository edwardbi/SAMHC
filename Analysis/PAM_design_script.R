# Compare PAM design with other designs
rm(list = ls())

folder_dir <- "./" # Change this to where the folder is
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
source("./MSAM1_slice.R")
source("./helper_funcs.R")

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

save_out_dir <- paste0(save_dir,"input_data_samples.RData")
saveRDS(list(cov = d_mat, ori_out = y_mat_null, gen_out = y_mat_alt), file = save_out_dir)

##################################################################
# PAM Result effect
##################################################################

# Hyperparameters
M_iter <- 5000
burnin <- 5000

p <- length(d_mat[1,1,])

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


post_result <- MCMC_SAM1(d_mat, mu0, sig0, a, b, kappa0, nu0, a1, b1, a2, b2, 
                         0.5, 0.3, 0.5, 10, M_iter, burnin)

save_PAM_MCMC_dir <- paste0(save_dir,"PAM_design_MCMC_application.RData")
saveRDS(post_result, file = save_PAM_MCMC_dir)

cls <- matrix(rep(NA, M_iter*sum(c(n1,n2_new,n3))), nrow = M_iter)
for(iter in (burnin+1):(burnin+M_iter)){
  z_mat_est <- post_result$all_z_mat[[iter]]
  z_mat_vec_est <- as.vector(t(z_mat_est))
  z_mat_vec_est <- z_mat_vec_est[!is.na(z_mat_vec_est)]
  cls[iter-burnin,] <- z_mat_vec_est
}
psm <- comp.psm(cls)
optmal_cluster <- minVI(psm, cls, method="all")

save_PAM_OptC_dir <- paste0(save_dir,"PAM_design_OptC_application.RData")
saveRDS(optmal_cluster, file = save_PAM_OptC_dir)

# Hyperparameters for PAM design
A <- n1 - n2_new
a0 <- 3
b0 <- 3
a0_y <- 3
b0_y <- 3
mu0_y <- 0
nu0_y <- 0.1

PAM_result_null <- PAM_design_pp_app(d_mat, y_mat_null, post_result, n1, n2_new, n3, 
                                     2, M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)

PAM_result_alt <- PAM_design_pp_app(d_mat, y_mat_alt, post_result, n1, n2_new, n3, 
                                    2, M_iter, burnin, A, a0, b0, mu0_y, nu0_y, a0_y, b0_y)

save_PAM_TE_null_dir <- paste0(save_dir,"PAM_design_TE_application_null.RData")
saveRDS(PAM_result_null, file = save_PAM_TE_null_dir)
save_PAM_TE_alt_dir <- paste0(save_dir,"PAM_design_TE_application_alt.RData")
saveRDS(PAM_result_alt, file = save_PAM_TE_alt_dir)

