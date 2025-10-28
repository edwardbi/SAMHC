library(cluster)
library(mvtnorm)
library(prodlim)
library(LaplacesDemon)

stirling <- function(nn){
  
  if (!exists("maxnn", where=parent.frame(2))){
    assign("maxnn", 1, envir=parent.frame(2))
    assign("allss", list(1), envir=parent.frame(2))
  }
  
  
  maxnn.local <- get("maxnn", envir=parent.frame(2))
  
  #only calculate this if needed
  if (nn > maxnn.local){
    allss.local <- get("allss", envir=parent.frame(2))
    
    allss.local[(length(allss.local) + 1):nn] <- 0
    
    for (mm in (maxnn.local + 1):nn){
      allss.local[[mm]] <- c(allss.local[[mm - 1]] * (mm - 1), 0) +
        c(0, allss.local[[mm - 1]])
      mss <- max(allss.local[[mm]])
      allss.local[[mm]] <- allss.local[[mm]] / mss
    }
    
    assign("maxnn", nn, envir=parent.frame(2))
    assign("allss", allss.local, envir=parent.frame(2))
  }
  
  assign("nn", nn, envir=parent.frame(2))
  ss <- eval(quote(allss[[nn]]), envir=parent.frame(2))
  return(ss)
}

# ConstWeight Function
const_weight <- function(w_p){
  
  m_wp <- 1 - w_p
  m_wp2 <- shift(cumprod(m_wp), 1)
  m_wp2[1] <- 1
  
  w <- w_p*m_wp2
  return(w)
  
}

xi_jk <- function(kappa,k){
  return((1-kappa)*kappa^(k-1))
}

stoc_trunc <- function(kappa, c_mat){
  c_mat_vec <- as.vector(t(c_mat))
  c_mat_vec <- c_mat_vec[!is.na(c_mat_vec)]
  
  u_ij <- runif(length(c_mat_vec), 0, xi_jk(kappa,c_mat_vec))
  K_c <- 1 + floor((log(u_ij) - log(1 - kappa))/log(kappa))
  J_c <- max(c_mat_vec)
  K_s <- max(c(K_c, J_c))
  xi_c <- xi_jk(kappa, 1:K_s)
  
  return(list(u_vec = u_ij, K_s = K_s, xi = xi_c))
  
}

# f_phi(x_ji, phi_k) computes prob of p(x_ji|phi_k)
f_phi <- function(x_ji, phi_k, sig_k){
  return(mvtnorm::dmvnorm(x_ji, mean = phi_k, sigma = sig_k))
}

# sample zi
sample_zji <- function(x_ji, u_ji, pi_j_p, xi_c, phi, sig){
  
  # Construct pi_j from pi_j_p
  pi_j <- const_weight(pi_j_p)
  
  #print(phi)
  
  k_prob <- rep(0, length(pi_j_p))
  for(k in 1:length(pi_j_p)){
    phi_k <- phi[k,]
    sig_k <- sig[k,,]
    pi_jk <- pi_j[k]
    if(pi_jk > 0 && (u_ji < xi_c[k])){
      k_prob[k] <- log(pi_jk) + log(f_phi(x_ji, phi_k, sig_k)) - log(xi_c[k])
    }else{
      k_prob[k] <- -Inf
    }
  }
  
  #print(k_prob)
  
  if(!is.finite(max(k_prob))){
    updated_zji <- sample(seq(1,length(pi_j_p)),1)
  }else{
    k_prob <- exp(k_prob - max(k_prob))
    updated_zji <- sample(seq(1,length(pi_j_p)),1,prob = k_prob)
  }
  
  return(updated_zji)
}

# Point to consider 1: is n_jk = 0 in p_j*?
# This model here is wrong.
sample_pi_jk_p <- function(jk_idx, z, alpha0, beta_p, p){
  
  # Construct beta
  beta <- const_weight(beta_p)
  
  # draw pi_jk_p
  z_j <- z[jk_idx[1],]
  if(jk_idx[2] %in% z_j){
    n_jk <- sum(z_j == jk_idx[2], na.rm = T)
  }else{
    n_jk <- 0
  }
  sum_n_js <- 0
  for(s in (jk_idx[2]+1):length(beta_p)){
    sum_n_js <- sum_n_js + sum(z_j == s, na.rm = T)
  }
  
  beta_dist_alpha <- alpha0*beta[jk_idx[2]]
  beta_dist_beta <- alpha0*(1-sum(beta[1:jk_idx[2]]))
  
  if(beta[jk_idx[2]] == 0){
    beta_dist_alpha <- alpha0*(1e-6)
  }
  
  if((1-sum(beta[1:jk_idx[2]])) == 0){
    beta_dist_beta <- alpha0*(1e-6)
  }
  
  # Revised version
  
  if(n_jk == 0){
    p_j_star_denum <- p[jk_idx[1]] + (1-p[jk_idx[1]])*
      (beta(beta_dist_alpha, beta_dist_beta)/beta(beta_dist_alpha, beta_dist_beta+sum_n_js))
    #print(c(p[jk_idx[1]], (1-p[jk_idx[1]]), beta(beta_dist_alpha, beta_dist_beta), 
    #        beta(beta_dist_alpha, beta_dist_beta+sum_n_js)))
    p_j_star <- p[jk_idx[1]]/p_j_star_denum
    if(runif(1) <= p_j_star){
      # temp solution for rbeta == 1
      rand_weight <- rbeta(1, beta_dist_alpha, beta_dist_beta+sum_n_js)
      if(rand_weight == 1){
        rand_weight <- rand_weight - 0.0001
      }
      return(rand_weight)
    }else{
      return(0)
    }
  }else{
    rand_weight <- rbeta(1, beta_dist_alpha+n_jk, beta_dist_beta+sum_n_js)
    if(rand_weight == 1){
      rand_weight <- rand_weight - 0.0001
    }
    return(rand_weight)
  }
  
}

sample_p_j <- function(j, pi_p, a, b){
  
  pi_j <- const_weight(pi_p[j,])
  num_pi_j_p_eq_zero <- sum(pi_p[j,] <= 0.001)
  
  p_j_upd <- rbeta(1, a + length(pi_j) - num_pi_j_p_eq_zero, b + num_pi_j_p_eq_zero)
  return(p_j_upd)
}

sample_beta_p <- function(pi_p, beta_p, alpha0, gamma, stepsize){
  
  const_beta <- function(beta_k_p, k, beta_p){
    beta_p_prop <- beta_p
    beta_p_prop[k] <- beta_k_p
    beta_prop <- const_weight(beta_p_prop)
    
    return(beta_prop)
  }
  
  log_likelihood_pi <- function(k, pi_p, beta, alpha0){
    all_ll_sum <- 0
    for(j in 1:length(pi_p[,1])){
      for(l in k:(length(pi_p[1,]))){
        if(pi_p[j,l] != 0){
          if(l == 1){
            all_ll_sum = all_ll_sum + log(dbeta(pi_p[j,1],alpha0*beta[1],alpha0*(1-beta[1])))
          }else{
            all_ll_sum = all_ll_sum + log(dbeta(pi_p[j,l],alpha0*beta[l],alpha0*(1-sum(beta[1:l]))))
          }
        }
      }
    }
    return(all_ll_sum)
  }
  
  log_prior <- function(beta_k_p, gamma){
    return(log(dbeta(beta_k_p,1,gamma)))
  }
  
  proposal <- function(beta_k_p_old, stepsize){
    beta_k_p_prop <- runif(1, min = beta_k_p_old - stepsize, max = beta_k_p_old + stepsize)
    if(beta_k_p_prop < 0){
      beta_k_p_prop <- -beta_k_p_prop
    }
    if(beta_k_p_prop > 1){
      beta_k_p_prop <- 2 - beta_k_p_prop
    }
    return(beta_k_p_prop)
  }
  
  MH_step <- function(k, beta_p, pi_p, alpha0, gamma, stepsize){
    beta_k_p <- beta_p[k]
    # proposal for beta_k_p
    proposal_betq_k_p <- proposal(beta_k_p, stepsize)
    
    # reconstruct beta
    beta_prop <- const_beta(proposal_betq_k_p, k, beta_p)
    beta_old <- const_beta(beta_k_p, k, beta_p)
    
    # acceptance ratio
    accept_ratio <- log_likelihood_pi(k, pi_p, beta_prop, alpha0) + log_prior(proposal_betq_k_p, gamma) - 
      log_likelihood_pi(k, pi_p, beta_old, alpha0) - log_prior(beta_k_p, gamma)
    
    ll_pi_p <- log_likelihood_pi(k, pi_p, beta_prop, alpha0)
    ll_pi_c <- log_likelihood_pi(k, pi_p, beta_old, alpha0)
    if(is.infinite(ll_pi_p) && is.infinite(ll_pi_c)){
      accept_ratio <- -Inf
    }
    
    if(log(runif(1)) <= accept_ratio){
      return(list(beta_k_p = proposal_betq_k_p, acc = 1))
    }else{
      return(list(beta_k_p = beta_k_p, acc = 0))
    }
  }
  
  updated_beta_p <- rep(NA, length(beta_p))
  accept_list <- rep(NA, length(beta_p))
  for(k in 1:(length(beta_p))){
    updated_beta_p_result <- MH_step(k, beta_p, pi_p, alpha0, gamma, stepsize)
    updated_beta_val <- updated_beta_p_result$beta_k_p
    if(updated_beta_val == 1){
      updated_beta_val <- updated_beta_val - 1e-4
    }
    updated_beta_p[k] <- updated_beta_val
    accept_list[k] <- updated_beta_p_result$acc
  }
  
  return(list(beta_p_upd = updated_beta_p, acc_list = accept_list))
  
}

sample_phi_k <- function(k, z_mat, x_mat, mu0, kappa0, nu0, sig0){
  
  dim_x_j <- rep(NA, (dim(x_mat)[1]))
  for(j in 1:(dim(x_mat)[1])){
    x_mat_p1_j <- x_mat[j,,1]
    if(!is.na(x_mat_p1_j[length(x_mat_p1_j)])){
      dim_x_j[j] <- length(x_mat_p1_j)
    }else{
      dim_x_j[j] <- which(is.na(x_mat_p1_j))[1] - 1
    }
  }
  
  item_dim <- x_mat[1,1,]
  cov_p <- dim(x_mat)[3] #length(item_dim)
  sum_x_ji <- rep(0, cov_p)
  n <- 0
  for(j in 1:(dim(x_mat)[1])){
    for(i in 1:dim_x_j[j]){#(dim(x_mat)[2])){
      if(z_mat[j,i] == k){
        sum_x_ji <- sum_x_ji + x_mat[j,i,]
        n <- n + 1
      }
    }
  }
  
  if(n > 0){
    avg_x_ji <- sum_x_ji/n
    sum_x_ji_sub_phi_sq <- matrix(rep(0, cov_p*cov_p), nrow = cov_p)
    for(j in 1:(dim(x_mat)[1])){
      for(i in 1:dim_x_j[j]){#(dim(x_mat)[2])){
        if(z_mat[j,i] == k){
          sum_x_ji_sub_phi_sq <- sum_x_ji_sub_phi_sq + ((x_mat[j,i,] - avg_x_ji) %*% t(x_mat[j,i,] - avg_x_ji))
        }
      }
    }
    
    kappa0_upd <- kappa0 + n
    nu0_upd <- nu0 + n
    mu0_upd <- (kappa0*mu0 + n*avg_x_ji)/kappa0_upd
    sig0_upd <- sig0 + sum_x_ji_sub_phi_sq + 
      ((kappa0*n)/kappa0_upd)*((avg_x_ji - mu0) %*% t(avg_x_ji - mu0))
    return(rnorminvwishart(n = 1, mu0_upd, kappa0_upd, sig0_upd, nu0_upd))
  }else{
    return(rnorminvwishart(n = 1, mu0, kappa0, sig0, nu0))
  }
}

sample_m_mat <- function(z_mat, beta_p, alpha0){
  
  beta <- const_weight(beta_p)
  
  m_mat <- matrix(rep(NA, length(z_mat[,1])*length(beta_p)), ncol = length(beta_p))
  
  for(j in 1:length(z_mat[,1])){
    z_j <- z_mat[j,!is.na(z_mat[j,])]
    #print(z_j)
    for(k in 1:length(beta_p)){
      n_jk <- sum(z_j == k)
      #print(k)
      #print(n_jk)
      if(n_jk > 1){
        prob_m <- rep(NA, n_jk)
        ss <- stirling(n_jk)
        for(m in 1:n_jk){
          prob_m[m] <- ss[m]*(alpha0*beta[k])^m
        }
        #print(prob_m)
        m_mat[j,k] <- sample(seq(1,n_jk), size = 1, prob = prob_m/sum(prob_m))
      }else if(n_jk == 1){
        m_mat[j,k] <- 1
      }else{
        m_mat[j,k] <- 0
      }
      #print(m_mat[j,k])
    }
  }
  
  return(m_mat)
  
}

sample_gamma <- function(m_mat, z_mat, gamma0, a, b){
  #z_mat, m_mat, gamma0, a, b){
  m_dd <- sum(as.vector(m_mat), na.rm = T)
  #K <- length(m_mat[1,])
  z_mat_vec <- as.vector(z_mat)
  z_mat_vec <- z_mat_vec[!is.na(z_mat_vec)]
  K <- length(unique(z_mat_vec))
  
  w <- rbeta(1, gamma0+1, m_dd)
  weights <- c((a+K-1)/(b-log(w)), m_dd)
  s <- rbinom(1, 1, prob = weights[1]/sum(weights))
  
  if(s == 1){
    gamma_upd <- rgamma(1, a+K, b-log(w))
  }else{
    gamma_upd <- rgamma(1, a+K-1, b-log(w))
  }
  
  return(gamma_upd)
}

sample_alpha0 <- function(z_mat, m_mat, alpha0, a, b){
  J <- length(z_mat[,1])
  w_vec <- rep(NA, J)
  s_vec <- rep(NA, J)
  for(j in 1:J){
    n_j <- length(z_mat[j,!is.na(z_mat[j,])])
    w_vec[j] <- rbeta(1, alpha0+1, n_j)
    s_vec[j] <- runif(1)*(alpha0+n_j) < n_j
  }
  m_dd <- sum(as.vector(m_mat), na.rm = T)
  alpha0_upd <- rgamma(1, a+m_dd-sum(s_vec), b - sum(log(w_vec)))
  return(alpha0_upd)
}

gibbs_SAM1 <- function(z_mat, pi_p_mat, p_vec, beta_p_vec, phi_vec, sig_mat, gamma_val, alpha0_val, 
                       x_mat, kappa, mu0, sig0, a, b, a1, b1, a2, b2, kappa0, nu0, stepsize, sd_alpha0){
  z_mat_upd <- z_mat
  pi_p_mat_upd <- pi_p_mat
  p_vec_upd <- p_vec
  beta_p_vec_upd <- beta_p_vec
  phi_vec_upd <- phi_vec
  sig_mat_upd <- sig_mat
  alpha0_upd <- alpha0_val
  gamma_upd <- gamma_val
  
  # varying dim of x_mat and z_mat
  dim_x_j <- rep(NA, length(z_mat_upd[,1]))
  for(j in 1:length(z_mat_upd[,1])){
    x_mat_p1_j <- x_mat[j,,1]
    if(!is.na(x_mat_p1_j[length(x_mat_p1_j)])){
      dim_x_j[j] <- length(x_mat_p1_j)
    }else{
      dim_x_j[j] <- which(is.na(x_mat_p1_j))[1] - 1
    }
  }
  
  # Stochastic Truncation and sample of U_ji
  stoc_trunc_iter <- stoc_trunc(kappa, z_mat_upd)
  u_vec <- stoc_trunc_iter$u_vec
  u_mat <- z_mat
  for(j in 1:length(z_mat_upd[,1])){
    if(j == 1){
      u_mat[1,1:dim_x_j[1]] <- u_vec[1:dim_x_j[1]]
    }else{
      u_mat[j,1:dim_x_j[j]] <- u_vec[(sum(dim_x_j[1:(j-1)])+1):sum(dim_x_j[1:j])]
    }
  }
  max_K <- stoc_trunc_iter$K_s
  xi_c <- stoc_trunc_iter$xi
  
  # Update all beta_k_p
  prev_K <- length(pi_p_mat_upd[1,])
  #print(phi_vec_upd)
  if(prev_K > max_K){
    pi_p_mat_upd <- pi_p_mat_upd[,1:max_K]
    beta_p_vec_upd <- beta_p_vec_upd[1:max_K]
    phi_vec_upd <- phi_vec_upd[1:max_K,]
    sig_mat_upd <- sig_mat_upd[1:max_K,,]
  }
  beta_p_result <- sample_beta_p(pi_p_mat_upd, beta_p_vec_upd, alpha0_upd, gamma_upd, stepsize)
  beta_p_vec_upd <- beta_p_result$beta_p_upd
  if(prev_K < max_K){
    beta_p_vec_upd <- c(beta_p_vec_upd, rbeta(max_K - prev_K, 1, gamma_upd))
    
    phi_vec_temp <- matrix(rep(NA, max_K*length(phi_vec_upd[1,])), ncol = length(phi_vec_upd[1,]))
    phi_vec_temp[1:prev_K,] <- phi_vec_upd
    sig_vec_temp <- array(dim = c(max_K, dim(sig_mat_upd)[2], dim(sig_mat_upd)[3]))
    sig_vec_temp[1:prev_K,,] <- sig_mat_upd
    for(idx_i in 1:(max_K-prev_K)){
      rand_NIW <- rnorminvwishart(n = 1, mu0, kappa0, sig0, nu0)
      phi_vec_temp[prev_K+idx_i,] <- rand_NIW$mu[1,]
      sig_vec_temp[prev_K+idx_i,,] <- rand_NIW$Sig
    }
    phi_vec_upd <- phi_vec_temp
    sig_mat_upd <- sig_vec_temp
    
    new_pi_p_mat <- matrix(rep(0, (max_K - prev_K)*length(pi_p_mat_upd[,1])), nrow = length(pi_p_mat_upd[,1]))
    pi_p_mat_upd <- cbind(pi_p_mat_upd, new_pi_p_mat)
  }
  
  #print(phi_vec_upd)
  
  #print(beta_p_vec_upd)
  #beta_p_vec_upd <- c(0.5, 1)
  #print("Done sample beta")
  
  #print(pi_p_mat_upd)
  #print(beta_p_vec_upd)
  
  # Update all pi_jk_p
  for(j in 1:length(pi_p_mat_upd[,1])){
    for(k in 1:max_K){
      pi_p_mat_upd[j,k] <- sample_pi_jk_p(c(j,k), z_mat_upd, alpha0_upd, beta_p_vec_upd, p_vec_upd)
    }
  }
  
  #print(pi_p_mat_upd)
  #pi_p_mat_upd <- matrix(c(0.5, 1, 0.999, 1, 0, 1), ncol = 2, byrow = TRUE)
  #print("Done sample pi_p_mat")
  
  # Update p
  for(j in 1:length(pi_p_mat_upd[,1])){
    p_vec_upd[j] <- sample_p_j(j, pi_p_mat_upd, a, b)
  }
  
  #print("Done sample p_vec")
  
  #print(phi_vec_upd)
  
  # Update all z_jis
  for(j in 1:length(z_mat_upd[,1])){
    pi_j_p <- pi_p_mat_upd[j,]
    #print(pi_j_p)
    for(i in 1:dim_x_j[j]){#length(z_mat_upd[1,])){
      x_ji <- x_mat[j,i,]
      u_ji <- u_mat[j,i]
      z_mat_upd[j,i] <- sample_zji(x_ji, u_ji, pi_j_p, xi_c, phi_vec_upd, sig_mat_upd)
    }
  }
  
  #print(z_mat_upd)
  
  #z_mat_upd <-c_mat
  
  #z_mat_upd <- z_mat_t
  
  #print("Done sample z_mat")
  
  # Update all phi_k
  #prev_phi_K <- length(phi_vec_upd[,1])
  #if(prev_phi_K > max_K){
  #  phi_vec_upd <- phi_vec_upd[1:max_K,]
  #  sig_mat_upd <- sig_mat_upd[1:max_K,,]
  #}
  #if(prev_phi_K < max_K){
  #  phi_vec_upd <- c(phi_vec_upd,rep(NA, max_K-prev_phi_K))
  #  sig_vec_upd <- c(sig_vec_upd,rep(NA, max_K-prev_phi_K))
  #}
  for(k in 1:max_K){
    sampled_point <- sample_phi_k(k, z_mat_upd, x_mat, mu0, kappa0, nu0, sig0)
    phi_vec_upd[k,] <- sampled_point$mu[1,]
    sig_mat_upd[k,,] <- sampled_point$Sigma
    #sample_phi_k(k, z_mat_upd, x_mat, sig_mat_upd, mu0, sig0)
  }
  
  m_mat <- sample_m_mat(z_mat_upd, beta_p_vec_upd, alpha0_upd)
  
  #print("Done sample m_mat")
  
  gamma_upd <- sample_gamma(m_mat, z_mat_upd, gamma_upd, a2, b2)
  
  #print("Done sample gamma")
  #print(gamma_upd)
  
  alpha0_upd <- sample_alpha0(z_mat_upd, m_mat, alpha0_upd, a1, b1)
  
  #print("Done sample alpha0")
  #print(alpha0_upd)
  
  return(list(z_mat_upd = z_mat_upd, pi_p_mat_upd = pi_p_mat_upd, 
              p_vec_upd = p_vec_upd, beta_p_vec_upd = beta_p_vec_upd, 
              beta_accpt = beta_p_result$acc_list, 
              phi_vec_upd = phi_vec_upd, sig_mat_upd = sig_mat_upd,
              gamma_upd = gamma_upd, alpha0_upd = alpha0_upd))
}

init_mats <- function(x_mat, K, a, b){
  
  #data_vec <- as.vector(t(x_mat))
  
  dim_x_j <- rep(NA, (dim(x_mat)[1]))#length(z_mat_upd[,1]))
  for(j in 1:(dim(x_mat)[1])){
    x_mat_p1_j <- x_mat[j,,1]
    if(!is.na(x_mat_p1_j[length(x_mat_p1_j)])){
      dim_x_j[j] <- length(x_mat_p1_j)
    }else{
      dim_x_j[j] <- which(is.na(x_mat_p1_j))[1] - 1
    }
  }
  
  
  # Assume data has all equal patients
  #memship_vec <- kmeans(data_vec, centers=K)
  #z_mat <- matrix(memship_vec$cluster, ncol = length(x_mat[1,]), byrow = TRUE)
  
  z_mat <- matrix(rep(NA, (dim(x_mat)[1])*(dim(x_mat)[2])), ncol = (dim(x_mat)[2]))
  
  for(j in 1:(dim(x_mat)[1])){
    z_mat[j,1:dim_x_j[j]] <- sample(seq(1,K), dim_x_j[j], replace = TRUE)
  }
  
  alpha0 <- 1
  gamma <- 1
  
  # Init phi_k
  phi_vec <- matrix(rep(NA, K*(dim(x_mat)[3])), ncol = (dim(x_mat)[3]))#rep(NA, K)
  sig_mat <- array(dim = c(K, dim(x_mat)[3], dim(x_mat)[3]))
  for(k in 1:K){
    sum_obs <- rep(0, dim(x_mat)[3])
    n_val <- 0
    #x_mat_obs <- NA
    for(j in 1:length(z_mat[,1])){
      for(i in 1:dim_x_j[j]){#length(z_mat[1,])){
        if(z_mat[j,i] == k){
          sum_obs <- sum_obs + x_mat[j,i,]
          n_val <- n_val + 1
          #if(is.na(x_mat_obs)){
          #  x_mat_obs <- x_mat[j,i,]
          #}else{
          #  x_mat_obs <- rbind(x_mat_obs, x_mat[j,i,])
          #}
        }
      }
    }
    phi_vec[k,] <- sum_obs/n_val
    sig_mat[k,,] <- diag(nrow = dim(x_mat)[3])
  }
  
  # Init beta_k'
  beta_p_vec <- rbeta(K, 1, gamma)
  
  
  # Init pi_jk'
  pi_p_mat <- matrix(rep(NA, length(z_mat[,1])*K), ncol = K)
  for(j in 1:length(z_mat[,1])){
    pi_p_mat[j,] <- rbeta(K, 1, alpha0)
  }
  
  # Init p_j
  p_vec <- rep(NA, length(z_mat[,1]))
  for(j in 1:length(z_mat[,1])){
    p_vec[j] <- rbeta(1, a, b)
  }
  
  return(list(z_mat = z_mat, pi_p_mat = pi_p_mat, beta_p_vec = beta_p_vec, 
              phi_vec = phi_vec, p_vec = p_vec, sig_mat = sig_mat,
              alpha0 = alpha0, gamma = gamma))
}

MCMC_SAM1 <- function(x_mat, mu0, sig0, a, b, kappa0, nu0, a1, b1, a2, b2, kappa,
                      stepsize, sd_alpha0, max_K, N_iter, burn_in){
  
  #(x_mat, n_clust, K, alpha0, gamma, method = "unif")
  init_params <- init_mats(x_mat, max_K, a, b)
  
  #print(init_params)
  
  z_mat_list <- list()
  pi_p_mat_list <- list()
  p_vec_list <- list()
  beta_p_vec_list <- list()
  phi_vec_list <- list()
  sig_mat_list <- list()
  gamma_list <- list()
  alpha0_list <- list()
  acc_list <- list()
  
  z_mat_list[[1]] <- init_params$z_mat
  pi_p_mat_list[[1]] <- init_params$pi_p_mat
  p_vec_list[[1]] <- init_params$p_vec
  beta_p_vec_list[[1]] <- init_params$beta_p_vec
  phi_vec_list[[1]] <- init_params$phi_vec
  sig_mat_list[[1]] <- init_params$sig_mat
  gamma_list[[1]] <- init_params$gamma
  alpha0_list[[1]] <- init_params$alpha0
  
  #print(pi_p_mat_list[[1]])
  
  for(iter in 1:(burn_in + N_iter)){
    if(iter %% 1000 == 0){
      print(iter)
    }
    
    #print(iter)
    
    gibbs_result <- gibbs_SAM1(z_mat_list[[iter]], pi_p_mat_list[[iter]], p_vec_list[[iter]], 
                               beta_p_vec_list[[iter]], phi_vec_list[[iter]], sig_mat_list[[iter]],
                               gamma_list[[iter]], alpha0_list[[iter]], x_mat, kappa, mu0, sig0,
                               a, b, a1, b1, a2, b2, kappa0, nu0, stepsize, sd_alpha0)
    z_mat_list[[iter+1]] <- gibbs_result$z_mat_upd
    pi_p_mat_list[[iter+1]] <- gibbs_result$pi_p_mat_upd
    p_vec_list[[iter+1]] <- gibbs_result$p_vec_upd
    beta_p_vec_list[[iter+1]] <- gibbs_result$beta_p_vec_upd
    phi_vec_list[[iter+1]] <- gibbs_result$phi_vec_upd
    sig_mat_list[[iter+1]] <- gibbs_result$sig_mat_upd
    gamma_list[[iter+1]] <- gibbs_result$gamma_upd
    alpha0_list[[iter+1]] <- gibbs_result$alpha0_upd
    acc_list[[iter]] <- gibbs_result$beta_accpt
  }
  
  return(list(all_z_mat = z_mat_list, all_pi_p_mat = pi_p_mat_list, 
              all_p_vec = p_vec_list, all_beta_vec = beta_p_vec_list,
              all_phi_vec = phi_vec_list, all_sig_mat = sig_mat_list,
              all_gamma = gamma_list, all_alpha0 = alpha0_list,
              all_acc = acc_list))
  
}
