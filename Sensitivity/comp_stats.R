library(mclust)
library(mcclust)
library(mcclust.ext)

comp_psm_optclust <- function(num_p_j, burnin, M_iter, post_z_mat, mat_flag = T, all_flag = T){
  
  cls <- matrix(rep(NA, M_iter*sum(num_p_j)), nrow = M_iter)
  for(iter in (burnin+1):(burnin+M_iter)){
    z_mat_est <- post_z_mat[[iter]]
    if(mat_flag == T){
      z_mat_vec_est <- as.vector(t(z_mat_est))
      z_mat_vec_est <- z_mat_vec_est[!is.na(z_mat_vec_est)]
    }else{
      z_mat_vec_est <- z_mat_est
    }
    cls[iter-burnin,] <- z_mat_vec_est
  }
  psm <- comp.psm(cls)
  if(all_flag == F){
    optmal_cluster <- minVI(psm, cls) 
  }else{
    optmal_cluster <- minVI(psm, cls, method="all", include.greedy = TRUE)
  }
  
  return(list(psm = psm, optclust = optmal_cluster))
  
}

nfd <- function(mat_a, mat_b){
  p <- length(mat_a[1,])
  sum_val <- 0
  for(i in 1:p){
    for(j in 1:p){
      sum_val <- sum_val + ((mat_a[i,j] - mat_b[i,j])^2)/(p^2)
    }
  }
  return(sum_val)
}

comp_randidx_ndf <- function(num_p_j, c_mat, psm, optmal_cluster, all_flag = T){
  
  c_mat_vec <- as.vector(t(c_mat))
  c_mat_vec <- c_mat_vec[!is.na(c_mat_vec)]
  if(all_flag == T){
    rand_idx <- adjustedRandIndex(optmal_cluster$cl[1,], c_mat_vec)
  }else{
    rand_idx <- adjustedRandIndex(optmal_cluster$cl, c_mat_vec)
  }
  
  all_mat_true <- matrix(rep(NA, (sum(num_p_j)*sum(num_p_j))), ncol = sum(num_p_j))
  c_mat_all_vec <- as.vector(t(c_mat))
  c_mat_all_vec <- c_mat_all_vec[!is.na(c_mat_all_vec)]
  for(i in 1:sum(num_p_j)){
    for(i1 in 1:sum(num_p_j)){
      if(c_mat_all_vec[i] == c_mat_all_vec[i1]){
        all_mat_true[i,i1] <- 1
      }else{
        all_mat_true[i,i1] <- 0
      }
    }
  }
  ndf_val <- nfd(psm, all_mat_true)
  
  return(list(randidx = rand_idx, ndf = ndf_val))
  
}
