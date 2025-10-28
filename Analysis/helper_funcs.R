# ConstWeight Function
const_weight <- function(w_p){
  
  m_wp <- 1 - w_p
  m_wp2 <- shift(cumprod(m_wp), 1)
  m_wp2[1] <- 1
  
  w <- w_p*m_wp2
  return(w)
  
}

rnorminvgamma <- function(n, mu, nu, alpha, beta){
  Y <- rinvgamma(n, alpha, beta)
  if(n == 1){
    X <- rnorm(1, mu, sd = sqrt(Y/nu))
    return(data.frame(sig = sqrt(Y), mu = X))
  }else{
    X <- rep(NA, n)
    for(i in 1:n){
      X[i] <- rnorm(1, mu, sd = sqrt(Y[i]/nu))
    }
    return(data.frame(sig = sqrt(Y), mu = X))
  }
}

# Posterior continuous outcome trt
post_trt_cont <- function(n, y1_vec, z1_vec, k, mu0, nu0, a0, b0){
  
  y_k_vec <- y1_vec[z1_vec == k]
  n1k <- length(y_k_vec)
  mu_post <- (nu0*mu0 + sum(y_k_vec))/(nu0 + n1k)
  nu_post <- nu0 + n1k
  a_post <- a0 + n1k/2
  b_post <- b0 + 0.5*(sum((y_k_vec - mean(y_k_vec))^2)) + 0.5*(n1k*nu0/(n1k+nu0))*(mean(y_k_vec) - mu0)^2
  return(rnorminvgamma(n = n, mu_post, nu_post, a_post, b_post))
  
}

# Posterior binary outcome trt
post_trt_binary <- function(n, y1_vec, z1_vec, k, a0, b0){
  
  y_k_vec <- y1_vec[z1_vec == k]
  sum_y <- sum(y_k_vec)
  n1k <- length(y_k_vec)
  a_p <- a0 + sum_y
  b_p <- b0 + n1k - sum_y
  return(rbeta(n,a_p,b_p))
  
}

# Posterior continuous outcome ctl
post_ctl_cont <- function(n, y2_vec, z2_vec, y3_vec, z3_vec, k, mu0, nu0, a0, b0, gamk){
  
  y2_k_vec <- y2_vec[z2_vec == k]
  sum_y2 <- sum(y2_k_vec)
  n2k <- length(y2_k_vec)
  y3_k_vec <- y3_vec[z3_vec == k]
  sum_y3 <- sum(y3_k_vec)
  n3k <- length(y3_k_vec)
  mu_pp <- (nu0*mu0 + gamk*sum(y3_k_vec))/(nu0 + gamk*n3k)
  nu_pp <- nu0 + gamk*n3k
  a_pp <- a0 + gamk*n3k/2
  b_pp <- b0 + (gamk/2)*(sum((y3_k_vec - mean(y3_k_vec))^2)) + 0.5*(gamk*n3k*nu0/(gamk*n3k+nu0))*(mean(y3_k_vec) - mu0)^2
  
  mu_post <- (nu_pp*mu_pp + sum(y2_k_vec))/(nu_pp + n2k)
  nu_post <- nu_pp + n2k
  a_post <- a_pp + n2k/2
  b_post <- b_pp + 0.5*(sum((y2_k_vec - mean(y2_k_vec))^2)) + 0.5*(n2k*nu_pp/(n2k+nu_pp))*(mean(y2_k_vec) - mu_pp)^2
  return(rnorminvgamma(n = n, mu_post, nu_post, a_post, b_post))
  
}

# Posterior binary outcome ctl
post_ctl_binary <- function(n, y2_vec, z2_vec, y3_vec, z3_vec, k, a0, b0, gamk){
  
  y2_k_vec <- y2_vec[z2_vec == k]
  sum_y2 <- sum(y2_k_vec)
  n2k <- length(y2_k_vec)
  y3_k_vec <- y3_vec[z3_vec == k]
  sum_y3 <- sum(y3_k_vec)
  n3k <- length(y3_k_vec)
  a_p <- a0 + sum_y2 + gamk*sum_y3
  b_p <- b0 + (n2k - sum_y2) + gamk*(n3k - sum_y3)
  return(rbeta(n,a_p,b_p))
  
}

comp_power_prior <- function(n, ps_result, ps_bor_rct, strata_list, mu0_y, nu0_y, a0, b0){
  
  # power prior
  trt_eff_cont <- 0
  trt_eff_cont_sig2 <- 0
  trt_eff_bin <- 0
  for(s in 1:5){
    gamk <- ps_bor_rct$Borrow[ps_bor_rct$Borrow[,1] == strata_list[s],11]
    # Continuous outcome
    theta_default_trt_cont <- post_trt_cont(n, ps_result$outcome_cont[ps_result$arm == 1], 
                                            ps_result$`_strata_`[ps_result$arm == 1], strata_list[s], 
                                            mu0_y, nu0_y, a0, b0)
    theta_default_trt_cont_mu <- theta_default_trt_cont$mu
    theta_default_trt_cont_sig <- (theta_default_trt_cont$sig)^2
    group3_out <- ps_result$outcome_cont[ps_result$group == 0 & ps_result$arm == 0]
    group3_strata <- ps_result$`_strata_`[ps_result$group == 0 & ps_result$arm == 0]
    group3_NA <- which(is.na(group3_strata))
    if(length(group3_NA) > 0){
      group3_out <- group3_out[-group3_NA]
      group3_strata <- group3_strata[-group3_NA]
    }
    theta_default_ctl_cont <- post_ctl_cont(n, ps_result$outcome_cont[ps_result$group == 1 & ps_result$arm == 0], 
                                            ps_result$`_strata_`[ps_result$group == 1 & ps_result$arm == 0], 
                                            group3_out, group3_strata, strata_list[s], mu0_y, nu0_y, a0, b0, gamk)
    theta_default_ctl_cont_mu <- theta_default_ctl_cont$mu
    theta_default_ctl_cont_sig <- (theta_default_ctl_cont$sig)^2
    trt_eff_cont <- trt_eff_cont + (0.2)*(theta_default_trt_cont_mu - theta_default_ctl_cont_mu)
    trt_eff_cont_sig2 <- trt_eff_cont_sig2 + (0.2^2)*(theta_default_trt_cont_sig + theta_default_ctl_cont_sig)
    # Binary outcome
    theta_default_trt_bin <- post_trt_binary(n, ps_result$outcome_bi[ps_result$arm == 1], 
                                             ps_result$`_strata_`[ps_result$arm == 1], strata_list[s], a0, b0)
    group3_out_bi <- ps_result$outcome_bi[ps_result$group == 0 & ps_result$arm == 0]
    if(length(group3_NA) > 0){
      group3_out_bi <- group3_out_bi[-group3_NA]
    }
    theta_default_ctl_bin <- post_ctl_binary(n, ps_result$outcome_bi[ps_result$group == 1 & ps_result$arm == 0], 
                                             ps_result$`_strata_`[ps_result$group == 1 & ps_result$arm == 0], 
                                             group3_out_bi, group3_strata, strata_list[s], a0, b0, gamk)
    #trt_eff_bin <- trt_eff_bin + (0.2)*(logit(theta_default_trt_bin) - logit(theta_default_ctl_bin))
    trt_eff_bin <- trt_eff_bin + (0.2)*(theta_default_trt_bin - theta_default_ctl_bin)
  }
  
  return(list(cont = trt_eff_cont, cont_sig2 = trt_eff_cont_sig2, bin = trt_eff_bin))
  
}

# Compute cluster-specific weighting factor
discount_factor <- function(w1, w2, z_vec, k, A, r = 2){
  w1_k <- w1[k]
  w2_k <- w2[k]
  #if(w1_k < (1/r)*w2_k){
  #  w1_k <- w2_k
  #}
  Ak <- round((r/(r-1))*(w1_k - (1/r)*w2_k)*A)
  nk <- sum(z_vec == k)
  gamma_k <- min(c(Ak,nk))/nk
  return(gamma_k)
  
}

# Merge singleton clusters
merge_unique <- function(k, common_clusters, all_phi, all_pi, z_vec){
  
  phi_k <- all_phi[k,]
  best_com_k <- NA
  best_dist <- NA
  for(kp in common_clusters){
    phi_kp <- all_phi[kp,]
    dist <- norm(phi_k, phi_kp)
    if(is.na(best_dist)){
      best_dist <- dist
      best_com_k <- kp
    }else if(dist < best_dist){
      best_dist <- dist
      best_com_k <- kp
    }
  }
  z_vec_upd <- z_vec
  z_vec_upd[z_vec_upd == k] <- best_com_k
  all_pi_upd <- all_pi
  all_pi_upd[best_com_k] <- all_pi_upd[best_com_k] + all_pi_upd[k]
  all_pi_upd[k] <- 0
  
  return(list(z_vec = z_vec_upd, pi_vec = all_pi_upd))
}

PAM_design_pp <- function(d_mat, y_cont_mat, y_bin_mat, post_result, n1, n2, n3, 
                          r, M_iter, burnin, A, a0, b0, mu0, nu0, a0n, b0n){
  
  # Hyperparameters
  A <- A
  a0 <- a0
  b0 <- b0
  mu0_y <- mu0
  nu0_y <- nu0
  a0_y <- a0n
  b0_y <- b0n
  
  # Overall treatment effect (Most interested in!)
  trt_eff_bin_MCMC <- rep(NA, M_iter)
  trt_eff_cont_MCMC <- rep(NA, M_iter)
  
  # Cluster-specific weights and treatment effect
  K <- 50
  pi_1_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  pi_2_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_trt_cont_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_aug_ctl_cont_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_trt_bin_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_aug_ctl_bin_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  
  # cluster-membership matrix
  cls <- matrix(rep(NA, M_iter*sum(c(n1,n2,n3))), nrow = M_iter)
  
  # In each MCMC of PAM:
  for(m in (burnin+1):(burnin+M_iter)){
    
    # The posterior cluster matrix from MCMC
    z_mat_est <- post_result$all_z_mat[[m]]
    g1_z <- z_mat_est[1,1:n1]
    g1_uniq_z <- unique(g1_z)
    g2_z <- z_mat_est[2,1:n2]
    g2_uniq_z <- unique(g2_z)
    g12_unique_z <- unique(c(g1_uniq_z,g2_uniq_z))
    g3_z <- z_mat_est[3,1:n3]
    g3_uniq_z <- unique(g3_z)
    
    # Save the cluster membership matrix as a row vector
    # 1:n1_actual = G1; (n1_actual+1):(n1+n2) = G2; (n1+n2+1):(n1+n2+n3) = G3
    z_mat_vec <- as.vector(t(z_mat_est))
    cls[(m-burnin),] <- z_mat_vec[!is.na(z_mat_vec)]
    
    # Step 2a) find common clusters of RWD with RCT
    common_clusts_G3 <- intersect(g12_unique_z, g3_uniq_z)
    d_mat_g3_common <- d_mat[3, g3_z %in% common_clusts_G3,]
    y_cont_vec_g3_common <- y_cont_mat[3, g3_z %in% common_clusts_G3]
    y_bin_vec_g3_common <- y_bin_mat[3, g3_z %in% common_clusts_G3]
    g3_z_common <- g3_z[g3_z  %in% common_clusts_G3]
    
    # Merge clusters of singleton
    common_clusts_G12 <- intersect(g1_uniq_z, g2_uniq_z)
    pi_p_mat <- post_result$all_pi_p_mat[[m]]
    all_pi_1 <- const_weight(pi_p_mat[1,])
    all_pi_2 <- const_weight(pi_p_mat[2,])
    all_phi <- post_result$all_phi_vec[[m]]
    
    # Merge unique clusters in Trt
    for(k1 in g1_uniq_z){
      if(!(k1 %in% common_clusts_G12)){
        merge_k1_result <- merge_unique(k1, common_clusts_G12, all_phi, all_pi_1, g1_z)
        all_pi_1 <- merge_k1_result$pi_vec
        g1_z <- merge_k1_result$z_vec
      }
    }
    
    # Merge unique clusters in Ctl
    for(k2 in g2_uniq_z){
      if(!(k2 %in% common_clusts_G12)){
        merge_k2_result <- merge_unique(k2, common_clusts_G12, all_phi, all_pi_2, g2_z)
        all_pi_2 <- merge_k2_result$pi_vec
        g2_z <- merge_k2_result$z_vec
      }
    }
    
    # Save weights for G1 and G2
    pi_1_MCMC[(m-burnin),1:length(pi_p_mat[1,])] <- all_pi_1
    pi_2_MCMC[(m-burnin),1:length(pi_p_mat[2,])] <- all_pi_2
    
    # Treatment effect
    trt_eff_bin <- 0
    trt_eff_cont <- 0
    for(k in common_clusts_G12){
      # Compute cluster k's treatment effect
      theta_bi_k_1 <- post_trt_binary(1, y_bin_mat[1,1:n1], g1_z, k, a0, b0)
      theta_cont_k_1_all <- post_trt_cont(1, y_cont_mat[1,1:n1], g1_z, k, mu0_y, nu0_y, a0_y, b0_y)
      theta_cont_k_1 <- theta_cont_k_1_all$mu
      
      if(k %in% common_clusts_G3){
        # When cluster k is common between G1, G2, and G3
        # Compute cluster-specific weight and weighting factor
        gamk <- discount_factor(all_pi_1, all_pi_2, g3_z_common, k, A, r = r)
        
        # Compute cluster k's control effect
        theta_bi_k_0 <- post_ctl_binary(1, y_bin_mat[2,1:n2], g2_z, 
                                        y_bin_vec_g3_common, g3_z_common, k, a0, b0, gamk)
        theta_cont_k_0_all <- post_ctl_cont(1, y_cont_mat[2,1:n2], g2_z, y_cont_vec_g3_common, 
                                            g3_z_common, k, mu0_y, nu0_y, a0_y, b0_y, gamk)
        theta_cont_k_0 <- theta_cont_k_0_all$mu
        
      }else{
        # When k is in G1 and G2, not in G3, no borrow
        theta_bi_k_0 <- post_trt_binary(1, y_bin_mat[2,1:n2], g2_z, k, a0, b0)
        theta_cont_k_0_all <- post_trt_cont(1, y_cont_mat[2,1:n2], g2_z, k, mu0_y, nu0_y, a0_y, b0_y)
        theta_cont_k_0 <- theta_cont_k_0_all$mu
      }
      
      # Overall treatment effect at cluster k
      
      trt_eff_bin <- trt_eff_bin + (theta_bi_k_1 - theta_bi_k_0)*all_pi_1[k] #(logit(theta_bi_k_1) - logit(theta_bi_k_0))*all_pi_1[k]
      trt_eff_cont <- trt_eff_cont + (theta_cont_k_1 - theta_cont_k_0)*all_pi_1[k]
      
      # Save cluster-specific trt eff
      theta_trt_cont_MCMC[(m-burnin), k] <- theta_cont_k_1
      theta_trt_bin_MCMC[(m-burnin), k] <- logit(theta_bi_k_1) #theta_bi_k_1
      theta_aug_ctl_cont_MCMC[(m-burnin), k] <- theta_cont_k_0
      theta_aug_ctl_bin_MCMC[(m-burnin), k] <- logit(theta_bi_k_0)#theta_bi_k_0
    }
    
    # Save overall treatment effect
    trt_eff_bin_MCMC[(m-burnin)] <- trt_eff_bin
    trt_eff_cont_MCMC[(m-burnin)] <- trt_eff_cont
    
    #if(m %% 100 == 0){
    #  print(trt_eff_cont)
    #}
  }
  
  return(list(p1 = pi_1_MCMC, p2 = pi_2_MCMC, 
              theta_trt_cont_k = theta_trt_cont_MCMC, theta_trt_bin_k = theta_trt_bin_MCMC, 
              theta_ctl_cont_k = theta_aug_ctl_cont_MCMC, theta_ctl_bin_k = theta_aug_ctl_bin_MCMC, 
              trt_eff_cont = trt_eff_cont_MCMC, trt_eff_bin = trt_eff_bin_MCMC))
  
}

PAM_design_pp_app <- function(d_mat, y_bin_mat, post_result, n1, n2, n3, r, M_iter, burnin, A, a0, b0, mu0, nu0, a0n, b0n){
  
  # Hyperparameters
  A <- A
  a0 <- a0
  b0 <- b0
  mu0_y <- mu0
  nu0_y <- nu0
  a0_y <- a0n
  b0_y <- b0n
  
  # Overall treatment effect (Most interested in!)
  trt_eff_bin_MCMC <- rep(NA, M_iter)
  
  # Cluster-specific weights and treatment effect
  K <- 50
  pi_1_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  pi_2_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_trt_bin_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_aug_ctl_bin_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  
  # cluster-membership matrix
  cls <- matrix(rep(NA, M_iter*sum(c(n1,n2,n3))), nrow = M_iter)
  
  # In each MCMC of PAM:
  for(m in (burnin+1):(burnin+M_iter)){
    
    # The posterior cluster matrix from MCMC
    z_mat_est <- post_result$all_z_mat[[m]]
    g1_z <- z_mat_est[1,1:n1]
    g1_uniq_z <- unique(g1_z)
    g2_z <- z_mat_est[2,1:n2]
    g2_uniq_z <- unique(g2_z)
    g12_unique_z <- unique(c(g1_uniq_z,g2_uniq_z))
    g3_z <- z_mat_est[3,1:n3]
    g3_uniq_z <- unique(g3_z)
    
    # Save the cluster membership matrix as a row vector
    # 1:n1_actual = G1; (n1_actual+1):(n1+n2) = G2; (n1+n2+1):(n1+n2+n3) = G3
    z_mat_vec <- as.vector(t(z_mat_est))
    cls[(m-burnin),] <- z_mat_vec[!is.na(z_mat_vec)]
    
    # Step 2a) find common clusters of RWD with RCT
    common_clusts_G3 <- intersect(g12_unique_z, g3_uniq_z)
    d_mat_g3_common <- d_mat[3, g3_z %in% common_clusts_G3,]
    y_bin_vec_g3_common <- y_bin_mat[3, which(g3_z %in% common_clusts_G3)]
    g3_z_common <- g3_z[g3_z  %in% common_clusts_G3]
    
    # Merge clusters of singleton
    common_clusts_G12 <- intersect(g1_uniq_z, g2_uniq_z)
    pi_p_mat <- post_result$all_pi_p_mat[[m]]
    all_pi_1 <- const_weight(pi_p_mat[1,])
    all_pi_2 <- const_weight(pi_p_mat[2,])
    all_phi <- post_result$all_phi_vec[[m]]
    
    # Merge unique clusters in Trt
    for(k1 in g1_uniq_z){
      if(!(k1 %in% common_clusts_G12)){
        merge_k1_result <- merge_unique(k1, common_clusts_G12, all_phi, all_pi_1, g1_z)
        all_pi_1 <- merge_k1_result$pi_vec
        g1_z <- merge_k1_result$z_vec
      }
    }
    
    # Merge unique clusters in Ctl
    for(k2 in g2_uniq_z){
      if(!(k2 %in% common_clusts_G12)){
        merge_k2_result <- merge_unique(k2, common_clusts_G12, all_phi, all_pi_2, g2_z)
        all_pi_2 <- merge_k2_result$pi_vec
        g2_z <- merge_k2_result$z_vec
      }
    }
    
    # Save weights for G1 and G2
    pi_1_MCMC[(m-burnin),1:length(pi_p_mat[1,])] <- all_pi_1
    pi_2_MCMC[(m-burnin),1:length(pi_p_mat[2,])] <- all_pi_2
    
    # Treatment effect
    trt_eff_bin <- 0
    for(k in common_clusts_G12){
      # Compute cluster k's treatment effect
      theta_bi_k_1 <- post_trt_binary(1, y_bin_mat[1,1:n1], g1_z, k, a0, b0)
      
      if(k %in% common_clusts_G3){
        # When cluster k is common between G1, G2, and G3
        # Compute cluster-specific weight and weighting factor
        gamk <- discount_factor(all_pi_1, all_pi_2, g3_z_common, k, A, r = r)
        
        # Compute cluster k's control effect
        theta_bi_k_0 <- post_ctl_binary(1, y_bin_mat[2,1:n2], g2_z, 
                                        y_bin_vec_g3_common, g3_z_common, k, a0, b0, gamk)
        
      }else{
        # When k is in G1 and G2, not in G3, no borrow
        theta_bi_k_0 <- post_trt_binary(1, y_bin_mat[2,1:n2], g2_z, k, a0, b0)
      }
      
      # Overall treatment effect at cluster k
      trt_eff_bin <- trt_eff_bin + (theta_bi_k_1 - theta_bi_k_0)*all_pi_1[k]
      
      # Save cluster-specific trt eff
      theta_trt_bin_MCMC[(m-burnin), k] <- theta_bi_k_1 #theta_bi_k_1
      theta_aug_ctl_bin_MCMC[(m-burnin), k] <- theta_bi_k_0 #theta_bi_k_0
    }
    
    # Save overall treatment effect
    trt_eff_bin_MCMC[(m-burnin)] <- trt_eff_bin
    
    #if(m %% 100 == 0){
    #  print(trt_eff_cont)
    #}
  }
  
  return(list(p1 = pi_1_MCMC, p2 = pi_2_MCMC, theta_trt_bin_k = theta_trt_bin_MCMC, 
              theta_ctl_bin_k = theta_aug_ctl_bin_MCMC, trt_eff_bin = trt_eff_bin_MCMC))
  
}

