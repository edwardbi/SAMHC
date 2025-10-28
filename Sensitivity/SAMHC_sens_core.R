PAM_design_pp_y_sens <- function(d_mat, y_cont_mat, post_result, n1, n2, n3, r,
                                 M_iter, burnin, A, a0, b0, mu0, nu0, a0n, b0n){
  
  # Hyperparameters
  A <- A
  a0 <- a0
  b0 <- b0
  mu0_y <- mu0
  nu0_y <- nu0
  a0_y <- a0n
  b0_y <- b0n
  
  # Overall treatment effect (Most interested in!)
  trt_eff_cont_MCMC <- rep(NA, M_iter)
  
  # Cluster-specific weights and treatment effect
  K <- 50
  pi_1_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  pi_2_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_trt_cont_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  theta_aug_ctl_cont_MCMC <- matrix(rep(NA, M_iter*K), ncol = K)
  
  # cluster-membership matrix
  cls <- matrix(rep(NA, M_iter*sum(c(n1,n2,n3))), nrow = M_iter)
  
  # In each MCMC of PAM:
  for(m in (burnin+1):(burnin+M_iter)){
    
    if(m %% 500 == 0){
      print(m)
    }
    
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
    trt_eff_cont <- 0
    for(k in common_clusts_G12){
      # Compute cluster k's treatment effect
      theta_cont_k_1_all <- post_trt_cont(1, y_cont_mat[1,1:n1], g1_z, k, mu0_y, nu0_y, a0_y, b0_y)
      theta_cont_k_1 <- theta_cont_k_1_all$mu
      
      if(k %in% common_clusts_G3){
        # When cluster k is common between G1, G2, and G3
        # Compute cluster-specific weight and weighting factor
        gamk <- discount_factor(all_pi_1, all_pi_2, g3_z_common, k, A, r = r)
        
        # Change to new idea
        y_g2_k_cont <- y_cont_mat[2, which(g2_z == k)]
        y_g3_k_cont <- y_cont_mat[3, which(g3_z == k)]
        #print(c(length(y_g2_k_cont), length(y_g3_k_cont)))
        if(length(y_g2_k_cont) > 2 && length(y_g3_k_cont) > 2){
          #print(c(length(y_g2_k_cont), length(y_g3_k_cont)))
          r_k_cont <- overlap(list(X1 = y_g2_k_cont, X2 = y_g3_k_cont), type = "1", pairsOverlap = F)
          gamk_cont <- min(gamk, r_k_cont$OV)
        }else{
          gamk_cont <- gamk
        }
        
        theta_cont_k_0_all <- post_ctl_cont(1, y_cont_mat[2,1:n2], g2_z, y_cont_vec_g3_common, 
                                            g3_z_common, k, mu0_y, nu0_y, a0_y, b0_y, gamk_cont)
        theta_cont_k_0 <- theta_cont_k_0_all$mu
        
      }else{
        # When k is in G1 and G2, not in G3, no borrow
        theta_cont_k_0_all <- post_trt_cont(1, y_cont_mat[2,1:n2], g2_z, k, mu0_y, nu0_y, a0_y, b0_y)
        theta_cont_k_0 <- theta_cont_k_0_all$mu
      }
      
      # Overall treatment effect at cluster k
      trt_eff_cont <- trt_eff_cont + (theta_cont_k_1 - theta_cont_k_0)*all_pi_1[k]
      
      # Save cluster-specific trt eff
      theta_trt_cont_MCMC[(m-burnin), k] <- theta_cont_k_1
      theta_aug_ctl_cont_MCMC[(m-burnin), k] <- theta_cont_k_0
    }
    
    # Save overall treatment effect
    trt_eff_cont_MCMC[(m-burnin)] <- trt_eff_cont
    
    #if(m %% 100 == 0){
    #  print(trt_eff_cont)
    #}
  }
  
  return(list(p1 = pi_1_MCMC, p2 = pi_2_MCMC, 
              theta_trt_cont_k = theta_trt_cont_MCMC,  
              theta_ctl_cont_k = theta_aug_ctl_cont_MCMC,  
              trt_eff_cont = trt_eff_cont_MCMC))
  
}

SAM_TTP <- function(d_mat, out_data, opt_clust, n1, n2, n3, r, a0, b0, a0_y, b0_y, mu0_y, nu0_y,
                    n_eff_sim = 10000){
  
  # Get the clusters
  g1_z <- opt_clust[1:n1]
  g1_uniq_z <- unique(g1_z)
  g2_z <- opt_clust[(n1+1):(n1+n2)]
  g2_uniq_z <- unique(g2_z)
  g12_unique_z <- unique(c(g1_uniq_z,g2_uniq_z))
  g3_z <- opt_clust[(n1+n2+1):(n1+n2+n3)]
  g3_uniq_z <- unique(g3_z)
  
  # Find common clusters of RWD and RCT
  common_clusts_G3 <- intersect(g12_unique_z, g3_uniq_z)
  d_mat_g3_common <- d_mat[3, g3_z %in% common_clusts_G3,]
  y_cont_vec_g3_common <- out_data[3, g3_z %in% common_clusts_G3]
  y_cont_vec_g2_common <- out_data[2, g2_z %in% common_clusts_G3]
  g3_z_common <- g3_z[g3_z  %in% common_clusts_G3]
  common_clusts_G12 <- intersect(g1_uniq_z, g2_uniq_z)
  
  trt_eff <- rep(0,n_eff_sim)
  for(k in common_clusts_G12){
    y_g1_k <- out_data[1, which(g1_z == k)]
    y_g2_k <- out_data[2, which(g2_z == k)]
    y_g3_k <- out_data[3, which(g3_z == k)]
    g3_d <- data.frame(y = mean(y_g3_k), y.se = sd(y_g3_k)/sqrt(length(y_g3_k)), 
                       n = length(y_g3_k), study = "Group 3")
    
    w1_k <- length(y_g1_k)/length(out_data[1, which(g1_z %in% common_clusts_G12)])
    w2_k <- length(y_g2_k)/length(out_data[1, which(g1_z %in% common_clusts_G12)])
    A <- n1-n2
    Ak <- round((r/(r-1))*(w1_k - (1/r)*w2_k)*A)
    nk <- sum(g3_z_common == k)
    gamk <- min(c(Ak,nk))/nk
    
    test_same <- ks.test(x = y_g2_k, y = y_g3_k)
    p_val <- test_same$p.value
    
    theta_cont_k_1_all <- post_trt_cont(n_eff_sim, out_data[1, 1:n1], 
                                        g1_z, k, mu0_y, nu0_y, a0_y, b0_y)
    theta_cont_k_1 <- theta_cont_k_1_all$mu
    
    if(k %in% common_clusts_G3 && p_val > 0.05){
      theta_cont_k_0_all <- post_ctl_cont(n_eff_sim, out_data[2,1:n2], 
                                          g2_z, y_cont_vec_g3_common, 
                                          g3_z_common, k, mu0_y, nu0_y, a0_y, b0_y, gamk)
      theta_cont_k_0 <- theta_cont_k_0_all$mu
    }else{
      theta_cont_k_0_all <- post_trt_cont(n_eff_sim, out_data[2,1:n2], 
                                          g2_z, k, mu0_y, nu0_y, a0_y, b0_y)
      theta_cont_k_0 <- theta_cont_k_0_all$mu
    }
    
    prop_k <- length(y_g1_k)/length(out_data[1, which(g1_z %in% common_clusts_G12)])
    trt_eff <- trt_eff + (theta_cont_k_1 - theta_cont_k_0)*prop_k
  }
  
  bias <- trt_eff - rep(3,n_eff_sim)
  mse <- mean(bias^2)
  
  return(list(trt_eff_cont = trt_eff, bias = bias, mse = mse))
  
}
