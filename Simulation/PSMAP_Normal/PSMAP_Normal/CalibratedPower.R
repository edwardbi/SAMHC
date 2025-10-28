###Functions ----
test_results = function(results, cutoff) {
  lower = summary(results[[1]], quantiles = cutoff)$quantiles
  upper = summary(results[[1]], quantiles = 1-cutoff)$quantiles
  reject = !((lower < 0 ) & (upper > 0))
  
  return(reject)
}

calib.power.type1 = function(list1, cutoff) {
  
  cutoff_results = sapply(list1,FUN = test_results, cutoff = cutoff)
  type1 = rowSums(cutoff_results)/dim(cutoff_results)[2]
  return(type1)
}

calib.power = function(results.type1, results.power, cutoff) {
  # Finds the calibrated power for each of the methods so type 1 error is fixed at approximately 5%
  
  # results.type1: MCMC chains when there is no treatment effect
  # results.power: MCMC chains when there is a treatment effect
  # cutoff: grid values to search over to find when type 1 error is roughly 5%
  
  
  results.return = list()
  
  all.data = list(cur = list(), pooled = list(), MAP = list(),
                  rMAP = list(), base = list(), NP = list(), pow = list(), mpow = list())
  
  for(i in 1:length(results.type1)) {
    all.data$cur[[i]] = results.type1[[i]]$Current
    all.data$pooled[[i]] = results.type1[[i]]$Pooled
    all.data$MAP[[i]] = results.type1[[i]]$MAP
    all.data$rMAP[[i]] = results.type1[[i]]$rMAP
    all.data$base[[i]] = results.type1[[i]]$BaSe
    all.data$NP[[i]] = results.type1[[i]]$NP
    all.data$pow[[i]] = results.type1[[i]]$PowerPrior
    all.data$mpow[[i]] = results.type1[[i]]$mPowerPrior
  }
  
  all.trt = list(cur = list(), pooled = list(), MAP = list(),
                 rMAP = list(), base = list(), NP = list(), pow = list(), mpow = list())
  
  for(i in 1:length(results.power)) {
    all.trt$cur[[i]] = results.power[[i]]$Current
    all.trt$pooled[[i]] = results.power[[i]]$Pooled
    all.trt$MAP[[i]] = results.power[[i]]$MAP
    all.trt$rMAP[[i]] = results.power[[i]]$rMAP
    all.trt$base[[i]] = results.power[[i]]$BSD
    all.trt$NP[[i]] = results.power[[i]]$NP
    all.trt$pow[[i]] = results.power[[i]]$PowerPrior
    all.trt$mpow[[i]] = results.power[[i]]$mPowerPrior
  }
  
  
  
  
  all.type1 = sapply(all.data,FUN = calib.power.type1, cutoff = cutoff)
  
  cutoff_inds = apply(all.type1,MARGIN = 2, FUN = function(x) max(which(x <= 0.05)))
  
  results.return$exact = apply(all.type1,MARGIN = 2, FUN = function(x) sum(x == 0.05))
  
  cutoff_inds[is.infinite(cutoff_inds)] = 1
  cutoff.power = cutoff[cutoff_inds]
  
  all.power = sapply(all.trt,FUN = calib.power.type1, cutoff = cutoff.power)
  
  tmp = diag(all.power)
  names(tmp) = colnames(all.power)
  
  results.return$power = tmp
  results.return$alltype1 = all.type1
  results.return$allpower = all.power
  return(results.return)
}

### Obtain calibrated power ----


data.frame(power = apply(rejection_results_trt,MARGIN = 1,FUN = sum)/dim(rejection_results_trt)[2],
           type1 = apply(rejection_results_no,MARGIN = 1,FUN = sum)/dim(rejection_results_trt)[2], row.names = names(overall_results_trt[[1]]))

data.frame(power = apply(rejection_results_trt3,MARGIN = 1,FUN = sum)/dim(rejection_results_trt3)[2],
           type1 = apply(rejection_results_no3,MARGIN = 1,FUN = sum)/dim(rejection_results_trt3)[2], row.names = names(overall_results_trt3[[1]]))


cutoff = seq(0.01, 0.05, length.out = 400)
tic()
tmp = calib.power(results.type1 = overall_results_no3,results.power = overall_results_trt3, cutoff = cutoff)
toc()

tic()
tmp1 = calib.power(results.type1 = overall_results_no,results.power = overall_results_trt, cutoff = cutoff)
toc()



