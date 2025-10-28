
library(RBesT)
library(coda)
library(rjags)
library(invgamma)
library(tictoc)
library(abind)
library(boot)
library(rstan)
library(MASS)
library(ggplot2)

source("ps-power_functions.R")
source("ps-power_stan.R")




test_results = function(results, cutoff) {
  reject = rep(NA,length(results))
  for(k in 1:length(results)) {
    
    lower = summary(results[[k]][,'effect'], quantiles = cutoff)$quantiles
    upper = summary(results[[k]][,'effect'], quantiles = 1-cutoff)$quantiles
    reject[k] = !((lower < 0 ) & (upper > 0))
  }
  
  return(reject)
}



# gen_X = function(n.hist.total, n.cur, mu,Sigma,alpha){
#   p=dim(mu)[1]
#   X.hist=NULL
#   X.cur=NULL
#   stop=FALSE
#   num_cur=0
#   num_hist=0
#   while(stop==FALSE){
#     X_initial=mvrnorm(1, mu=mu, Sigma = Sigma)
#     X_i= X_convert(t(as.matrix(X_initial)), c(1, 2, 3, 4)) #convert the first four covariates to binary
#     ps_i = t(alpha) %*% X_i
#     prob_i=inv.logit(ps_i)
#     Z_i= rbern(1,prob_i)
#     if(Z_i==1 & num_cur<n.cur){
#       X.cur=rbind(X.cur,t(X_i))
#       num_cur=dim(X.cur)[1]
#     }else if(Z_i==0 & num_hist<n.hist.total){
#       X.hist=rbind(X.hist,t(X_i))
#       num_hist=dim(X.hist)[1]
#     }
#     
#     if(num_cur==n.cur&num_hist==n.hist.total){
#       stop=TRUE
#     }
#   }
#   return(list(X.hist,X.cur))
# }



X_convert = function(X, binary_col) {
  #Convert a particular set of continuous covariates to binary
  
  for (i in 1:length(binary_col)) {
    current_col = X[, i]
    X[which(current_col > 0), i] = 1
    X[which(current_col <= 0), i] = 0
  }
  return(t(X))
}


gen_Y = function(X, beta, var_eps,trt.effect=0,study.effect=0) {
  #Generate the outcome variable given the covariates
  N = dim(X)[1]
  Y = matrix(0, N)
  for (i in 1:N) {
    xb = t(beta) %*% X[i, ]
    Y[i, ] = xb +trt.effect+ study.effect+rnorm(1, 0, var_eps)
  }
  return(Y)
}




#summarize results
summarize_results = function(results, theta.true,cutoff,parameter) {
  n.methods=length(results)-1
  theta_mean = rep(NA,n.methods)
  theta_bias= rep(NA,n.methods)
  theta_MSE= rep(NA,n.methods)
  theta_CI_width= rep(NA,n.methods)
  theta_coverage= rep(NA,n.methods)
  reject= rep(NA,n.methods)
  for(k in 1:n.methods) {
    theta_draws=results[[k]][,parameter]
    theta_mean[k] = mean(theta_draws[[1]])
    theta_bias[k] = mean(theta_draws[[1]]) - theta.true
    theta_MSE[k] = mean((theta_draws[[1]] - theta.true) ^ 2)
    lower = summary(theta_draws, quantiles = cutoff)$quantiles
    upper = summary(theta_draws, quantiles = 1-cutoff)$quantiles
    
    
    theta_CI_width[k] = upper - lower
    theta_coverage[k] = ifelse(lower < theta.true &
                                 upper > theta.true, 1, 0)
    reject[k]= ifelse((lower < 0 ) & (upper > 0),0,1)
    
    
  }
  results=list(mean=theta_mean,bias=theta_bias,MSE=theta_MSE,CI_width=theta_CI_width,CI_coverage=theta_coverage,reject=reject)
  return(results)
}




MAP.prior = function(ybar.hist, SE.hist, prior.std, parameters = "theta.new", niter, sigma) {
  # Obtain prior samples and ESS when using the MAP prior, returns the MCMC object contain the parameters
  # listed in "parameters"
  
  # ybar.hist = mean from each of the historical trials
  # SE.hist = std error from each of the historical trials
  # prior.std = hyper-parameter, std dev of half normal placed on tau
  # parameters = name of parameters to return choice of: 
  ## mu: overall mean for random effects on control
  ## tau: the precision for random effects for control i.e. between trial precision
  ## theta: the individual random effects
  ## theta.new: a new sample from the prior 
  
  n.hist = length(ybar.hist)
  tau.hat = (1/ SE.hist^2)
  
  dataTemp = list(
    'ybar' = ybar.hist,
    'n.hist' = n.hist,
    'std_heterogeneity' = prior.std,
    'tau.hat' = tau.hat
  )
  
  model <-
    jags.model(
      file = "MAP_Prior.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  MAP_model <-
    coda.samples(model,
                 variable.names = parameters,
                 thin = 5,
                 n.iter = niter, progress.bar = 'none')
  
  
  mix.res <- automixfit(MAP_model[[1]][,1], Nc=1:4, thresh=0, type="norm")
  ess.res=ess( mix.res, method="elir",sigma=sigma)
  return(list(MAP_model,ess.res))
}




#' MAP.fit = function(ybar.trt, sigma.hat.trt, n.trt, ybar.control, sigma.hat.control, n.con, 
#'                    ybar.hist, SE.hist, prior.std , alpha.sigma = 1, beta.sigma = 1, parameters = "theta", niter) {
#'   # Obtain posterior samples when using the MAP prior, returns the MCMC object contain the parameters
#'   # listed in "parameters"
#'   
#'   # sigma.hat.trt/sigma.hat.control = MLE estimate of variance for trt/control
#'   # ybar.hist = mean from each of the historical trials
#'   # SE.hist = std error from each of the historical trials
#'   # prior.std = hyper-parameter, std dev of half normal placed on tau
#'   # alpha.sigma / beta.sigma = hyper-parameter, alpha/beta for InverseGamma placed on var on trt/control
#'   
#'   # parameters = name of parameters to return choice of: 
#'   ## effects: the difference between treatment mean and current control mean
#'   ## mu: overall mean for random effects on control
#'   ## mu.trt: mean for the treatment
#'   ## tau: the precision for random effects for control i.e. between trial precision
#'   ## tau.control/tau.trt: precision for control/treatment i.e. within trial precision
#'   ## theta: the individual random effects
#'   
#'   n.hist = length(ybar.hist)
#'   tau.hat = (1/ SE.hist^2)
#'   
#'   dataTemp = list(
#'     'ybar' = ybar.hist,
#'     'n.hist' = n.hist,
#'     'std_heterogeneity' = prior.std,
#'     'tau.hat' = tau.hat,
#'     #'ybar.trt' = ybar.trt,
#'     'ybar.control' = ybar.control,
#'     'sigma.hat.control' = sigma.hat.control,
#'     # 'sigma.hat.trt' = sigma.hat.trt,
#'     # 'n.trt' = n.trt,
#'     'n.con' = n.con,
#'     'alpha.sigma' =  alpha.sigma,
#'     'beta.sigma' = beta.sigma,
#'     'y1' = ybar.control
#'     # 'y2' = ybar.trt
#'   )
#'   
#'   model <-
#'     jags.model(
#'       file = "MAPNormal.txt",
#'       data = dataTemp,
#'       n.chains = 1,
#'       n.adapt = 0.2 * niter,
#'       quiet = TRUE
#'     )
#'   update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
#'   MAP_model <-
#'     coda.samples(model,
#'                  variable.names = parameters,
#'                  thin = 5,
#'                  n.iter = niter, progress.bar = 'none')
#'   
#'   return(list(MAP_model,ess.res))
#' }



MAP.fit2 = function(ybar.trt, sigma.hat.trt, n.trt, ybar.control, sigma.hat.control, n.con, 
                    ybar.hist, SE.hist, prior.std, alpha.sigma = 1, beta.sigma = 1, parameters = "effect", niter) {
  # Obtain posterior samples when using the MAP prior, returns the MCMC object contain the parameters
  # listed in "parameters"
  
  # sigma.hat.trt/sigma.hat.control = MLE estimate of variance for trt/control
  # ybar.hist = mean from each of the historical trials
  # SE.hist = std error from each of the historical trials
  # prior.std = hyper-parameter, std dev of half normal placed on tau
  # alpha.sigma / beta.sigma = hyper-parameter, alpha/beta for InverseGamma placed on var on trt/control
  
  # parameters = name of parameters to return choice of: 
  ## effects: the difference between treatment mean and current control mean
  ## mu: overall mean for random effects on control
  ## mu.trt: mean for the treatment
  ## tau: the precision for random effects for control i.e. between trial precision
  ## tau.control/tau.trt: precision for control/treatment i.e. within trial precision
  ## theta: the individual random effects
  
  n.hist = length(ybar.hist)
  tau.hat = (1/ SE.hist^2)
  
  dataTemp = list(
    'ybar' = ybar.hist,
    'n.hist' = n.hist,
    'std_heterogeneity' = prior.std,
    'tau.hat' = tau.hat,
    'ybar.trt' = ybar.trt,
    'ybar.control' = ybar.control,
    'sigma.hat.control' = sigma.hat.control,
    'sigma.hat.trt' = sigma.hat.trt,
    'n.trt' = n.trt,
    'n.con' = n.con,
    'alpha.sigma' =  alpha.sigma,
    'beta.sigma' = beta.sigma,
    'y1' = ybar.control,
    'y2' = ybar.trt
  )
  
  model <-
    jags.model(
      file = "MAPNormal2.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  MAP_model <-
    coda.samples(model,
                 variable.names = parameters,
                 thin = 5,
                 n.iter = niter, progress.bar = 'none')
  
  return(MAP_model)
}



Power.fit = function(ybar.trt, sigma.hat.trt, n.trt, ybar.control, sigma.hat.control, n.con, 
                     hist.ybar, SE.hist,n.hist.samples,prior.trt.var = 10^6, alpha.sigma = 1, beta.sigma = 1,
                     alpha = 0.1, niter, return.var = FALSE) {
  # Fits the power prior and returns an MCMC object containing the mean of treatment/control and difference,
  # can return the estimate of the treatment/control variance using 'return.var'
  
  # ybar.trt/sigma.trt/n.trt: the sample mean, MLE variance (over n not n-1) and sample size for treatement
  # ybar.control/sigma.hat.control/n.con: the sample mean, MLE variance (over n not n-1) and sample size for control
  # alpha.sigma/beta.sigma: parameters for the IG prior on treatment/control variance
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # hist.ybar/SE.hist/n.hist.samples: vector of historical means/standard errors/sample sizes
  # alpha: power parameter 
  # return.var: return the estimate of the observation-level variance for the treatment/control 
  
  ybar.hist = mean(hist.ybar)
  sum.sq = sum((SE.hist^2)*(n.hist.samples)*(n.hist.samples - 1) + n.hist.samples*(hist.ybar^2))
  sigma.hist = (sum.sq -  sum(n.hist.samples)*(ybar.hist^2))/(sum(n.hist.samples) - 1)
  sigma.hat.con = sigma.hat.control
  
  burn.in = (0.2)*niter
  # Create vectors to store parameters
  mu.trt = rep(NA,niter+burn.in)
  var.trt = rep(NA,niter+burn.in)
  mu = rep(NA,niter+burn.in)
  var.control = rep(NA,niter+burn.in)
  diff.effect = rep(NA,niter+burn.in)
  
  # Specify hyper-parameters
  prior.var = prior.trt.var 
  alpha.var = alpha.sigma
  beta.var = beta.sigma
  
  
  # Intialize the parameters
  mu.trt[1] = rnorm(1,sd = sqrt(prior.var))
  var.trt[1] = invgamma::rinvgamma(n=1,shape = alpha.var,rate = beta.var)
  mu[1] = rnorm(1,sd = sqrt(prior.var)) # Hyper-prior for control mean is actually jefferys prior i.e propto constant
  var.control[1] = invgamma::rinvgamma(n=1,shape = alpha.var,rate = beta.var)
  diff.effect[1] = mu.trt[1] - mu[1]
  
  for( i in 2:(niter+burn.in)) {
    # Update treatment mean
    mean.mu.trt = (n.trt * prior.var * ybar.trt)/(n.trt * prior.var + var.trt[i-1])
    var.mu.trt = (var.trt[i-1]*prior.var)/(n.trt * prior.var + var.trt[i-1])
    mu.trt[i] = rnorm(n = 1,mean = mean.mu.trt,sd = sqrt(var.mu.trt))
    
    # Update treatment variance
    a.trt = n.trt/2 + alpha.var
    b.trt = (n.trt/2)*(sigma.hat.trt + ((ybar.trt - mu.trt[i])^2)) + beta.var
    var.trt[i] =invgamma:: rinvgamma(n = 1,shape = a.trt,rate = b.trt)
    
    # Update control mean
    mean.mu.con = ((sigma.hist * n.con * ybar.control) + (var.control[i-1] * alpha * sum(n.hist.samples) * ybar.hist))/(n.con*sigma.hist + var.control[i-1] * alpha* sum(n.hist.samples))
    var.mu.con = (var.control[i-1] * sigma.hist)/(n.con * sigma.hist + var.control[i-1] * alpha * sum(n.hist.samples))
    mu[i] = rnorm(n = 1,mean = mean.mu.con,sd = sqrt(var.mu.con))
    
    # Update control variance
    a.con = n.con/2 + alpha.var
    b.con = (n.con/2)*(sigma.hat.con + ((ybar.control - mu[i])^2)) + beta.var
    var.control[i] = invgamma::rinvgamma(n = 1,shape = a.con,rate = b.con)
    
    diff.effect[i] = mu.trt[i] - mu[i]
  }
  
  # if(return.var) {
  #   Power.model = mcmc(data.frame(mu = mu[-(1:burn.in)],mu.trt = mu.trt[-(1:burn.in)],effect = diff.effect[-(1:burn.in)],var.trt =  var.trt[-(1:burn.in)],var.con = var.control[-(1:burn.in)]))
  #   
  # } else {
  #   
  #   Power.model = mcmc(data.frame(mu = mu[-(1:burn.in)],mu.trt = mu.trt[-(1:burn.in)],effect = diff.effect[-(1:burn.in)]))
  # }
  
  diff.effect= diff.effect[-(1:burn.in)]
  diff.effect= diff.effect[seq(1,length(diff.effect),by=5)]
  Power.model=mcmc(data.frame(effect = diff.effect))
  return(mcmc.list(Power.model))
}











PS_MAP.fit = function(tau.init, target.ESS, sigma, n.cur,ybar.hist, SE.hist,overlap_coefficients, niter,lim, data.indiv) {
  # Obtain posterior MCMC samples when using the PS-MAP prior
  # target.ESS = pre-specified sample size
  # sigma = the fixed reference scale
  # n.cur = number of subjects in each stratum of the current trial
  # ybar.hist = mean from each of the historical trials /strata
  # SE.hist = std error from each of the historical trials /strata
  # overlap_coefficients = overlap coefficients from PS stratum
  # alpha.sigma / beta.sigma = hyper-parameter, alpha/beta for InverseGamma placed on var on trt/control
  
  # parameters = name of parameters to return choice of: 
  ## effects: the difference between treatment mean and current control mean
  ## mu: overall mean for random effects on control
  ## mu.trt: mean for the treatment
  ## tau: the precision for random effects for control i.e. between trial precision
  ## tau.control/tau.trt: precision for control/treatment i.e. within trial precision
  ## theta: the individual random effects
  ##adjust.coefficients=(median(overlap_coefficients)/overlap_coefficients)*mean(SE.hist)
  n.strata = length(ybar.hist)
  total.cur=sum(n.cur)
  WS1=n.cur/total.cur
  ess.res=c(0)
  tau.res=c(0)
  direction=0
  # ess.res.strata=NULL
  tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.init #hyperparameter for HN(.)
  stop=FALSE
  low=lim[1]
  high=lim[2]
  tau.scale=tau.init
  
  #while(stop==FALSE){
    theta.pred=list()
    # ess.res.str=rep(0,n.strata)
    
    
  #  if((high-low<=0.01)){
  #    stop=TRUE
  #    tau.scale=tau.res[which.min(abs(ess.res-target.ESS))]
      # print(ess.res)
      # print(tau.res)
      # print(which.min(abs(ess.res-target.ESS)))
      # print(tau.scale)
   #   tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.scale
      
   # }
    
    
    for(i in 1:n.strata){
      
      dataTemp = list(
        'ybar' = ybar.hist[i],
        'n.hist' = 1,
        'std_heterogeneity' = tau.prior[i],
        'tau.hat' = (1/ SE.hist[i]^2)
      )
      
      model <-
        jags.model(
          file = "MAP_Prior.txt",
          data = dataTemp,
          n.chains = 1,
          n.adapt = 0.2 * niter,
          quiet = TRUE
        )
      update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
      MAP_model <-
        coda.samples(model,
                     variable.names = "theta.new",
                     thin = 5,
                     n.iter = niter, progress.bar = 'none')
      
      theta.pred[[i]]=c(MAP_model[[1]][,1])
      
      #mix.res <- automixfit(theta.pred[[i]], Nc=1:4, thresh=0, type="norm")
      
      #ess.res.str[i]=ess( mix.res , method="elir",sigma=sigma)
      
    }
    
    
    theta.pred = do.call(cbind, theta.pred)
    theta=WS1%*%t(theta.pred)
    
    mix.res <- automixfit(theta[1,], Nc=1:4, thresh=0, type="norm")
    
    ess.res=c(ess.res,ess(mix.res , method="elir",sigma=sigma))
    tau.res=c(tau.res,tau.scale)
    
    
    
  #   if(abs(ess.res[length(ess.res)]-target.ESS)<5){
  #     direction=0
  #     stop=TRUE
  #   }else if(ess.res[length(ess.res)]>target.ESS){
  #     direction=-1
  #     low=tau.scale
  #     tau.scale=(low+high)/2
  #     tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.scale
  #   }else{
  #     direction=1
  #     high=tau.scale
  #     tau.scale=(low+high)/2
  #     tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.scale
  #   }
  #   
  # }
  
  
  #print(tau.res[length(tau.res)])
  posterior.indiv<-postmix(mix.res, data.indiv)
  results=list(mix.res,ess.res[length(ess.res)],posterior.indiv,tau.res[length(tau.res)])
  return(results)
}






PS_MAP.fit2=function(ybar.trt, sigma.hat.trt, n.trt, tau.init, target.ESS, sigma, n.cur,ybar.hist, SE.hist,overlap_coefficients,
                     lim, data.indiv,prior.trt.var = 10^6, alpha.sigma = 1, beta.sigma = 1,niter) {
  # Fits the PS-MAP prior and returns an MCMC object containing the mean of treatment/control and difference
  
  # ybar.trt/sigma.trt/n.trt: the sample mean, MLE variance (over n not n-1) and sample size for treatement
  # ybar.control/sigma.hat.control/n.con: the sample mean, MLE variance (over n not n-1) and sample size for control
  # alpha.sigma/beta.sigma: parameters for the IG prior on treatment/control variance
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # hist.ybar/SE.hist/n.hist.samples: vector of historical means/standard errors/sample sizes
  
  
  burn.in = (0.2)*niter
  # Create vectors to store parameters
  mu.trt = rep(NA,niter+burn.in)
  var.trt = rep(NA,niter+burn.in)
  mu = rep(NA,niter+burn.in)
  var.control = rep(NA,niter+burn.in)
  diff.effect = rep(NA,niter+burn.in)
  
  # Specify hyper-parameters
  prior.var = prior.trt.var 
  alpha.var = alpha.sigma
  beta.var = beta.sigma
  
  
  
  # Intialize the parameters
  mu.trt[1] = rnorm(1,sd = sqrt(prior.var))
  var.trt[1] = rinvgamma(n=1,shape = alpha.var,rate = beta.var)
  mu[1] = rnorm(1,sd = sqrt(prior.var)) # Hyper-prior for control mean is actually jefferys prior i.e propto constant
  diff.effect[1] = mu.trt[1] - mu[1]
  
  for( i in 2:(niter+burn.in)) {
    # Update treatment mean
    mean.mu.trt = (n.trt * prior.var * ybar.trt)/(n.trt * prior.var + var.trt[i-1])
    var.mu.trt = (var.trt[i-1]*prior.var)/(n.trt * prior.var + var.trt[i-1])
    mu.trt[i] = rnorm(n = 1,mean = mean.mu.trt,sd = sqrt(var.mu.trt))
    
    # Update treatment variance
    a.trt = n.trt/2 + alpha.var
    b.trt = (n.trt/2)*(sigma.hat.trt + ((ybar.trt - mu.trt[i])^2)) + beta.var
    var.trt[i] = rinvgamma(n = 1,shape = a.trt,rate = b.trt)
    
  }
  
  PS_MAP_fit= PS_MAP.fit(tau.init,target.ESS,sigma, n.cur, ybar.hist, SE.hist, overlap_coefficients , 
                         niter=100000,lim,data.indiv)
  
  PS_MAP_ESS=PS_MAP_fit[[2]] 
  PS_MAP_post=PS_MAP_fit[[3]]
  tau_scale_limits=PS_MAP_fit[[4]]
  
  mu.con=rmix(PS_MAP_post,n=niter)
  mu.con= mu.con[seq(1,length(mu.con),by=5)]
  mu.trt=mu.trt[-(1:burn.in)]
  mu.trt= mu.trt[seq(1,length(mu.trt),by=5)]
  diff.effect = mu.trt - mu.con[1:length(mu.trt)]
  
  #PS_MAP.model = mcmc(data.frame(mu = mu.con[1:length(mu.trt)],mu.trt = mu.trt,effect = diff.effect))
  PS_MAP.model = mcmc(data.frame(effect = diff.effect))
  
  fit.result=list(mcmc.list(PS_MAP.model),PS_MAP_ESS,tau_scale_limits)
  return(fit.result)
}




print("loading stan model...")

model <- stan_model(model_code = powerps, model_name = "powerps")

print("done")

PS_Power.fit=function(ps_data,outcomes,S=5,A){
  
  
  ana.ovl  <- rwePSDist(ps_data)
  
  rS <- ana.ovl$Dist[1:S] / sum(ana.ovl$Dist[1:S])
  
  Y = outcomes
  dta = cbind(Y, ps_data$data)
  
  theta_draws <-
    rwePsPowDrawPost(
      dta,
      v.outcome = "Y",
      A = A,
      RS = rS,
      type = "continuous"
    )$post.theta
  
  
  
  results=theta_draws
  return(results)
}







PS_Power.fit2 = function(ybar.trt, sigma.hat.trt, n.trt, ps_data,outcomes,A,prior.trt.var = 10^6, alpha.sigma = 1, beta.sigma = 1,
                         niter) {
  # Fits the PS-power prior and returns an MCMC object containing the mean of treatment/control and difference,
  # can return the estimate of the treatment/control variance using 'return.var'
  
  # ybar.trt/sigma.trt/n.trt: the sample mean, MLE variance (over n not n-1) and sample size for treatement
  # ybar.control/sigma.hat.control/n.con: the sample mean, MLE variance (over n not n-1) and sample size for control
  # alpha.sigma/beta.sigma: parameters for the IG prior on treatment/control variance
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # hist.ybar/SE.hist/n.hist.samples: vector of historical means/standard errors/sample sizes
  # alpha: power parameter 
  # return.var: return the estimate of the observation-level variance for the treatment/control 
  
  
  burn.in = (0.2)*niter
  # Create vectors to store parameters
  mu.trt = rep(NA,niter+burn.in)
  var.trt = rep(NA,niter+burn.in)
  mu = rep(NA,niter+burn.in)
  var.control = rep(NA,niter+burn.in)
  diff.effect = rep(NA,niter+burn.in)
  
  # Specify hyper-parameters
  prior.var = prior.trt.var 
  alpha.var = alpha.sigma
  beta.var = beta.sigma
  
  
  
  # Intialize the parameters
  mu.trt[1] = rnorm(1,sd = sqrt(prior.var))
  var.trt[1] = rinvgamma(n=1,shape = alpha.var,rate = beta.var)
  mu[1] = rnorm(1,sd = sqrt(prior.var)) # Hyper-prior for control mean is actually jefferys prior i.e propto constant
  
  
  for( i in 2:(niter+burn.in)) {
    # Update treatment mean
    mean.mu.trt = (n.trt * prior.var * ybar.trt)/(n.trt * prior.var + var.trt[i-1])
    var.mu.trt = (var.trt[i-1]*prior.var)/(n.trt * prior.var + var.trt[i-1])
    mu.trt[i] = rnorm(n = 1,mean = mean.mu.trt,sd = sqrt(var.mu.trt))
    
    # Update treatment variance
    a.trt = n.trt/2 + alpha.var
    b.trt = (n.trt/2)*(sigma.hat.trt + ((ybar.trt - mu.trt[i])^2)) + beta.var
    var.trt[i] = rinvgamma(n = 1,shape = a.trt,rate = b.trt)
    
    
  }
  
  mu.con=PS_Power.fit(ps_data,outcomes,S=5,A)
  mu.con= mu.con[seq(1,length(mu.con),by=5)]
  mu.trt=mu.trt[-(1:burn.in)]
  mu.trt= mu.trt[seq(1,length(mu.trt),by=5)]
  diff.effect = mu.trt - mu.con[1:length(mu.trt)]
  
  #Power.model = mcmc(data.frame(mu = mu.con[1:length(mu.trt)],mu.trt = mu.trt,effect = diff.effect))
  
  Power.model = mcmc(data.frame(effect = diff.effect))
  return(mcmc.list(Power.model))
}



PS_regroup<-function(X.hist.combine,Y.hist.combine,X.cur.con,Y.cur.con,p){
  #Stratification based on propensity scores
  S=5
  Y.combine=rbind(Y.hist.combine,Y.cur.con)
  dta = rbind(X.hist.combine, X.cur.con)
  covs = NULL
  for (k in 1:p) {
    covs[k] = paste("V", k, sep = "")
  }
  colnames(dta) = covs
  group = c(rep(0, dim(X.hist.combine)[1]), rep(1, dim(X.cur.con)[1])) #historical data =0; current control=1
  dta = cbind(group, dta)
  dta = as.data.frame(dta)
  ## estimate PS by logistic regression
  ana.ps <-
    rwePS(
      dta,
      v.grp = "group",
      v.covs = covs,
      nstrata = S,
      type = "logistic"
    )
  
  ## get overlapping coefficients
  ana.ovl  <- rwePSDist(ana.ps)
  
  ## get rs
  rS <- ana.ovl$Dist[1:S] 
  
  
  
  X.hist.stratum =list()
  y.hist.stratum=list()
  n.hist.stratum=NULL
  X.cur.stratum =list()
  y.cur.stratum=list()
  n.cur.stratum=NULL
  for(i in 1:S) {
    index=which((ana.ps$data$group==0)&(ana.ps$data$"_strata_"==i))
    n.hist.stratum[i]=length(index)
    X.hist.stratum [[i]] <- ana.ps$data[index,2:11]
    y.hist.stratum[[i]]<-Y.combine[index]
    
    indey2=which((ana.ps$data$group==1)&(ana.ps$data$"_strata_"==i))
    n.cur.stratum[i]=length(indey2)
    X.cur.stratum [[i]] <- ana.ps$data[indey2,2:11]
    y.cur.stratum[[i]]<-Y.combine[indey2]
  }
  
  stratum.ybar = sapply(y.hist.stratum,FUN = mean)
  SE.stratum = sapply(y.hist.stratum, FUN = sd)/sqrt(n.hist.stratum)
  stratum.ybar.cur = sapply(y.cur.stratum,FUN = mean)
  SE.stratum.cur = sapply(y.cur.stratum, FUN = sd)/sqrt(n.cur.stratum)
  sigma.stratum.cur=sapply(y.cur.stratum, FUN = sd)*(n.cur.stratum-1)/n.cur.stratum
  out=list(rS=rS,X.hist=X.hist.stratum, Y.hist=y.hist.stratum, Ybar.hist=stratum.ybar,
           SE.hist= SE.stratum,Ybar.cur=stratum.ybar.cur,SE.cur=SE.stratum.cur,sigma.cur=sigma.stratum.cur,
           n.hist=n.hist.stratum, n.cur=n.cur.stratum,PS=ana.ps)
  return(out)
  
}




# 
# PS_regroup_v<-function(X.hist.combine,Y.hist.combine,X.cur.con,Y.cur.con, S=5,p){
#   Y.combine=as.matrix(c(Y.hist.combine,Y.cur.con))
#   dta = as.matrix(c(X.hist.combine, X.cur.con))
#   covs = NULL
#   for (k in 1:p) {
#     covs[k] = paste("V", k, sep = "")
#   }
#   colnames(dta) = covs
#   group = c(rep(0, length(X.hist.combine)), rep(1, length(X.cur.con))) #historical data =0; current control=1
#   dta = cbind(group, dta)
#   dta = as.data.frame(dta)
#   ## estimate PS by logistic regression
#   ana.ps <-
#     rwePS(
#       dta,
#       v.grp = "group",
#       v.covs = covs,
#       nstrata = S,
#       type = "logistic"
#     )
#   
#   ## get overlapping coefficients
#   ana.ovl  <- rwePSDist(ana.ps)
#   
#   ## get rs
#   rS <- ana.ovl$Dist[1:S] 
#   
#   
#   
#   X.hist.stratum =list()
#   y.hist.stratum=list()
#   n.hist.stratum=NULL
#   X.cur.stratum =list()
#   y.cur.stratum=list()
#   n.cur.stratum=NULL
#   for(i in 1:S) {
#     index=which((ana.ps$data$group==0)&(ana.ps$data$"_strata_"==i))
#     n.hist.stratum[i]=length(index)
#     X.hist.stratum [[i]] <- ana.ps$data[index,2]
#     y.hist.stratum[[i]]<-Y.combine[index]
#     
#     indey2=which((ana.ps$data$group==1)&(ana.ps$data$"_strata_"==i))
#     n.cur.stratum[i]=length(indey2)
#     X.cur.stratum [[i]] <- ana.ps$data[indey2,2]
#     y.cur.stratum[[i]]<-Y.combine[indey2]
#   }
#   
#   stratum.ybar = sapply(y.hist.stratum,FUN = mean)
#   SE.stratum = sapply(y.hist.stratum, FUN = sd)/sqrt(n.hist.stratum)
#   sigma=mean(sapply(y.hist.stratum, FUN = sd))
#   stratum.ybar.cur = sapply(y.cur.stratum,FUN = mean)
#   SE.stratum.cur = sapply(y.cur.stratum, FUN = sd)/sqrt(n.cur.stratum)
#   sigma.stratum.cur=sapply(y.cur.stratum, FUN = sd)*(n.cur.stratum-1)/n.cur.stratum
#   out=list(rS=rS,X.hist=X.hist.stratum, Y.hist=y.hist.stratum, Ybar.hist=stratum.ybar,
#            SE.hist= SE.stratum,Ybar.cur=stratum.ybar.cur,SE.cur=SE.stratum.cur,sigma.cur=sigma.stratum.cur,
#            n.hist=n.hist.stratum, n.cur=n.cur.stratum,PS=ana.ps,sigma=sigma)
#   return(out)
#   
# }















# Single_arm_simulation = function(trt.effect=0, niter, true.mu.hist=1.2, true.mu.cur=1, sd.mu.hist = sqrt(1.5),sd.mu.cur=1,MAP_ESS=90, Power_A=90) {
#   ### Generate Data ----
#   
#   
#   # Historic
#   n.hist = 1 # number of historical trials
#   n.hist.total=3000 # total number of patients 
#   n.hist.samples = rep(n.hist.total/n.hist,n.hist) # number of patients in each trial
#   
#   
#   # Current 
#   n.con = 200 #number of patients in control
#   n.trt = 200 #number of patients in treatment
#   
#   
#   # Generate covariates X
#   p = 10 #dimension of the covariates
#   rho = 0.1 #correlations of the covariates
#   
#   
#   beta = matrix(1, p)
#   
#   #mu.hist =cbind(rep(1.2,p),rep(1.3,p),rep(1.4,p),rep(1.5,p),rep(1.8,p)) #assume the mean of covariates are the different across historical trials
#   mu.hist = matrix(true.mu.hist, p)
#   mu.cur = matrix(true.mu.cur, p) # mean of covariates from the current study
#   S0 = diag(x = rep(sd.mu.hist^2, p))
#   S0[which(S0 == 0)] = rho * (sd.mu.hist^2)
#   S1 = diag(x = rep(sd.mu.cur^2, p))
#   S1[which(S1 == 0)] = rho * (sd.mu.cur^2)
#   
#   
#   
#   
#   hist.list = list()
#   X.hist =list()
#   for(i in 1:n.hist) {
#     
#     X.hist.init = mvrnorm( n.hist.samples[i], mu =mu.hist[,i], Sigma = S0)
#     X.hist[[i]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
#     
#     # }
#     # X.hist =list()
#     # for(i in 1:n.hist){
#     #   n=c(0,cumsum(n.hist.samples))
#     #   X.hist.combine= X.hist.combine[sample(nrow(X.hist.combine)),]
#     #   X.hist[[i]]=X.hist.combine[(n[i]+1):n[i+1],]
#     
#     hist.list[[i]] = gen_Y(X.hist[[i]], beta, var_eps=1)
#     
#   }
#   sigma=mean(sapply(hist.list, FUN = sd))
#   X.hist.combine=abind(X.hist, along=1)
#   Y.hist.combine=abind(hist.list, along=1)
#   hist.ybar = sapply(hist.list,FUN = mean)
#   SE.hist = sapply(hist.list, FUN = sd)/sqrt(n.hist.samples)
#   
#   
#   
#   # Current 
#   #generate covariates for control
#   X.cur.init <- mvrnorm( n.con, mu =mu.cur, Sigma = S1) #assume trt and control from the same population
#   X.cur.con = t(X_convert(X.cur.init, c(1, 2, 3, 4))) #convert the first four covariates to binary
#   
#   #X.cur.con=x.cur[1:n.con,] #covaraites of the control group
#   #X.cur.trt=x.cur[(n.con+1):(n.trt+n.con),] #covariates of the treatment group
#   
#   
#   #generate outcomes for control
#   y.control = gen_Y(X.cur.con, beta, var_eps=1)
#   ybar.control = mean(y.control)
#   sigma.hat.control = var( y.control) * (n.con - 1) / n.con
#   
#   
#   #generate covariates for treatment
#   X.cur.init.trt <- mvrnorm( n.trt, mu =mu.cur, Sigma = S1) #assume trt and control from the same popultion
#   x.cur.trt = t(X_convert(X.cur.init.trt, c(1, 2, 3, 4))) #convert the first four covariates to binary
#   
#   #generate outcomes for treatment
#   y.treatment = gen_Y(x.cur.trt, beta, var_eps=1,trt.effect)
#   ybar.trt = mean( y.treatment)
#   sigma.hat.trt = var( y.treatment) * (n.trt - 1) / n.trt
#   
#   
#   
#   #Calculate propensity scores
#   
#   regroup_data<-PS_regroup(X.hist.combine,Y.hist.combine,X.cur.con,y.control,p=p)
#   
#   
#   #print(regroup_data$rS)
#   #print(regroup_data$rS/sum(regroup_data$rS))
#   
#   ### Fit models ----
#   
#   # current_model = current.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#   #                             ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                             alpha.sigma = 1,beta.sigma = 1,parameters = 'effect',niter = niter)
#   # 
#   # pooled_model = pooled.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#   #                           ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                           hist.ybar = hist.ybar,SE.hist = SE.hist,n.hist.samples = n.hist.samples,
#   #                           alpha.sigma = 1,beta.sigma = 1,parameters = 'effect',niter = niter)
#   # 
#   # 
#   
#   # 
#   # PS_MAP_model = PS_MAP.fit(ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                           ybar.hist = regroup_data$Ybar.hist,SE.hist =  regroup_data$SE.hist,
#   #                           overlap_coefficients = regroup_data$rS,alpha.sigma = 1,beta.sigma = 1,prior.std = 5,parameters = "theta", niter = niter)
#   # 
#   
#   # PS_MAP_model_v1 = PS_MAP.fit.v1(ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                                 ybar.hist = regroup_data$Ybar.hist,SE.hist =  regroup_data$SE.hist,
#   #                                 overlap_coefficients = regroup_data$rS,alpha.sigma = 1,beta.sigma = 1,parameters = "theta", niter = niter)
#   # 
#   
#   # PS_MAP_model  = PS_MAP.fit.v2( regroup_data$Ybar.cur, regroup_data$SE.cur, regroup_data$n.cur,  ybar.hist = regroup_data$Ybar.hist, SE.hist =  regroup_data$SE.hist,
#   #                                   overlap_coefficients = regroup_data$rS,alpha.sigma = 1, beta.sigma = 1,prior.std=5, parameters = "theta", niter)
#   # 
#   # PS_MAP_model  = PS_MAP.fit.v3( regroup_data$Ybar.cur, regroup_data$SE.cur, regroup_data$n.cur,  ybar.hist = regroup_data$Ybar.hist, SE.hist =  regroup_data$SE.hist,
#   #                                overlap_coefficients = regroup_data$rS, parameters = "theta", niter)
#   # 
#   # PS_MAP_model_v1 = PS_MAP.fit.v4( regroup_data$Ybar.cur, regroup_data$SE.cur,regroup_data$sigma.cur, regroup_data$n.cur,  ybar.hist = regroup_data$Ybar.hist, SE.hist =  regroup_data$SE.hist,
#   #                  overlap_coefficients = regroup_data$rS,alpha.sigma = 1, beta.sigma = 1, parameters = "theta", niter)
#   # 
#   # 
#   # 
#   
#   PS_MAP_results  =PS_MAP.fit(tau.scale=1,A=MAP_ESS,sigma=sigma, n.cur=regroup_data$n.cur, ybar.hist = regroup_data$Ybar.hist, SE.hist =  regroup_data$SE.hist,overlap_coefficients = regroup_data$rS, 
#                               niter,lim=c(0.00001,2),data.indiv = y.control)
#   
#  
#   MAP_model = MAP.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#                       ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#                       ybar.hist = hist.ybar,SE.hist = SE.hist,
#                       prior.std = 5,alpha.sigma = 1,beta.sigma = 1,parameters = "theta", niter = niter)
#   
#   rMAP_model = rMAP.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#                         ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#                         ybar.hist = hist.ybar,SE.hist = SE.hist,w = 0.9,
#                         prior.std = 5,alpha.sigma = 1,beta.sigma = 1,parameters = "theta", niter = niter)
#   
#   PS_Power_results = PS_Power.fit(regroup_data$PS, rbind(Y.hist.combine,y.control), S=5,A=Power_A,sigma=sigma)
#   PS_Power_model=PS_Power_results[[1]]
#   # BaSe_model = BaSe.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#   #                       ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                       ybar.hist = hist.ybar,SE.hist = SE.hist,
#   #                       prior.prec = 1/5,alpha.sigma = 1,beta.sigma = 1,
#   #                       Nclus = 10,alpha.DPP = 1,parameters = "effect",niter = niter)
#   # 
#   # Power_model = Power.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#   #                         ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                         hist.ybar = hist.ybar,SE.hist = SE.hist,n.hist.samples = n.hist.samples,
#   #                         prior.trt.var = 10^5,alpha.sigma = 1,beta.sigma = 1,alpha = 0.1,
#   #                         niter = niter,return.var = FALSE)
#   # 
#   # mPower_model = mPower.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
#   #                           ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
#   #                           hist.ybar = hist.ybar,SE.hist = SE.hist,n.hist.samples = n.hist.samples,
#   #                           prior.trt.var = 10^5,alpha.sigma = 1,beta.sigma = 1,
#   #                           alpha.delta = 1,beta.delta = 1,niter = niter,return.var = FALSE)
#   # 
#   # Gather Results 
#   # results = list(Current = NULL, Pooled = NULL, MAP = NULL,
#   #               rMAP = NULL, BSD = NULL, PowerPrior = NULL, mPowerPrior = NULL)
#   
#   results=list(PSMAP=NULL,MAP=NULL,rMAP=NULL,PSPower=NULL,regroup_data=NULL,PSMAPess=NULL)
#   
#   #results$Current = current_model
#   #results$Pooled = pooled_model
#   
#   results$PSMAP=PS_MAP_model
#   
#   results$MAP = MAP_model
#   results$rMAP = rMAP_model
#   results$PSPower = PS_Power_model
#   results$regroup_data=regroup_data
#   results$PSMAPess=PS_MAP_results[[2]]
#   #results$PSPoweress=PS_Power_results[[2]]
#   # results$BSD = BaSe_model
#   # results$PowerPrior = Power_model
#   # results$mPowerPrior = mPower_model
#   
#   return(results)
#   
# }







Two_arm_simulation = function(trt.effect=0, niter, MAP_ESS=90, Power_A=90, prior.std=0.5, alpha=0.1,limits=c(0.001,2),getlimits=FALSE) {
  ### Generate Data ----
  
  
  # Historic
  n.hist = 3 # number of historical trials
  n.hist.total=600 # total number of patients 
  n.hist.samples = rep(n.hist.total/n.hist,n.hist) # number of patients in each trial
  
  
  # Current 
  n.trt = 300 #number of patients in treatment
  n.con = 150 #number of patients in control
  
  #Generate covariates X
  p = 10 #dimension of the covariates
  rho = 0.1 #correlations of the covariates
  eps.var=1
  beta_0 = rep(1,p)
  
  mu.hist = cbind(rep(1,p),rep(0.5,p),rep(1.2,p))
  S0 = diag(x = rep(1^2, p))
  S0[which(S0 == 0)] = rho * (1^2)
  
  
  mu.cur = matrix(1, p) # mean of covariates from the current study
  S1 = diag(x = rep(1^2, p))
  S1[which(S1 == 0)] = rho * (1^2)
  
  
  hist.list = list()
  X.hist =list()
  delta=c(0,0.2,-0.1)
  
  # for(i in 1:(n.hist-1)) {
  #   X.hist.init = mvrnorm(n.hist.samples[i], mu =mu.hist[,i], Sigma = S0)
  #   X.hist[[i]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
  #   hist.list[[i]] = gen_Y(X.hist[[i]], beta=beta_0, var_eps=eps.var,study.effect=delta[i])
  # }
  # 
  #   mix.index=sum(rbinom(n.hist.samples[n.hist],1,1/3))
  #   X.hist.init_1 = mvrnorm(mix.index, mu =mu.hist[,n.hist-1], Sigma = S0)
  #   X.hist.init_2 = mvrnorm(n.hist.samples[n.hist]-mix.index, mu =mu.hist[,n.hist-2], Sigma = S0)
  #   X.hist.init=rbind(X.hist.init_1,X.hist.init_2)
  #   X.hist[[n.hist]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
  #   
  
  
  for(i in 1:n.hist) {
    X.hist.init = mvrnorm(n.hist.samples[i], mu =mu.hist[,i], Sigma = S0)
    X.hist[[i]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
    hist.list[[i]] = gen_Y(X.hist[[i]], beta=beta_0, var_eps=eps.var,study.effect=delta[i])
  }
  
  
  sigma=mean(sapply(hist.list, FUN=sd))
  X.hist.combine=abind(X.hist, along=1)
  Y.hist.combine=abind(hist.list, along=1)
  hist.ybar = sapply(hist.list,FUN = mean)
  SE.hist = sapply(hist.list, FUN = sd)/sqrt(n.hist.samples)
  
  
  
  # Current 
  #generate covariates for current study
  X.cur.init = mvrnorm((n.con+n.trt), mu =mu.cur, Sigma = S1)#assume trt and control from the same population
  X.cur = t(X_convert(X.cur.init, c(1, 2, 3, 4)))
  X.cur.con=X.cur[1:n.con,] #covaraites of the control group
  X.cur.trt=X.cur[(n.con+1):(n.trt+n.con),] #covariates of the treatment group
  
  
  #generate outcomes for control
  y.control = gen_Y(X.cur.con, beta=beta_0, var_eps=eps.var)
  ybar.control = mean(y.control)
  sigma.hat.control = var( y.control) * (n.con - 1) / n.con
  
  
  #generate outcomes for treatment
  y.treatment = gen_Y(X.cur.trt, beta=beta_0, var_eps=eps.var,trt.effect=trt.effect)
  ybar.trt = mean( y.treatment)
  sigma.hat.trt = var( y.treatment) * (n.trt - 1) / n.trt
  
  
  
  #Calculate propensity scores
  
  regroup_data<-PS_regroup(X.hist.combine,Y.hist.combine,X.cur.con,y.control,p=p)
  
  
  
  
  ### Fit models ----
  
  PS_MAP_fit=PS_MAP.fit2(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt, n.trt=n.trt,
                         tau.init=prior.std,target.ESS=MAP_ESS,sigma=sigma,
                         n.cur=regroup_data$n.cur, ybar.hist = regroup_data$Ybar.hist, 
                         SE.hist =  regroup_data$SE.hist,overlap_coefficients = regroup_data$rS, 
                         lim=limits,data.indiv = y.control,
                         prior.trt.var = 10^6, alpha.sigma = 1, beta.sigma = 1, niter = niter)
  
  PS_MAP_model=PS_MAP_fit[[1]]
  PS_MAP_ESS= PS_MAP_fit[[2]]
  
  MAP_model = MAP.fit2(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
                       ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
                       ybar.hist = hist.ybar,SE.hist = SE.hist,
                       prior.std = 5, alpha.sigma = 1,beta.sigma = 1, parameters = "effect", niter = niter)
  
  MAP_ESS = MAP.prior(ybar.hist= hist.ybar,SE.hist = SE.hist, prior.std = prior.std, parameters = "theta.new", niter, sigma=sigma)[[2]]
  
  
  PS_Power_model = PS_Power.fit2(ybar.trt=ybar.trt, sigma.hat.trt = sigma.hat.trt, n.trt=n.trt, ps_data=regroup_data$PS,
                                 outcomes=c(Y.hist.combine,y.control),
                                 A=Power_A, prior.trt.var = 10^6, alpha.sigma = 1, beta.sigma = 1, niter=niter)
  
  Power_model = Power.fit(ybar.trt = ybar.trt,sigma.hat.trt = sigma.hat.trt,n.trt = n.trt,
                          ybar.control = ybar.control,sigma.hat.control = sigma.hat.control,n.con = n.con,
                          hist.ybar = hist.ybar,SE.hist = SE.hist,n.hist.samples = n.hist.samples,
                          prior.trt.var = 10^5,alpha.sigma = 1,beta.sigma = 1,alpha = alpha,
                          niter = niter,return.var = FALSE)
  
  
  # Gather Results 
  
  results=list(PSMAP=NULL,MAP=NULL,PSPower=NULL,PowerPrior=NULL,ESS=NULL)
  #results=list(PSMAP=NULL,MAP=NULL,ESS=NULL)
  
  
  results$PSMAP=PS_MAP_model
  results$MAP = MAP_model
  results$PSPower = PS_Power_model
  results$PowerPrior = Power_model
  results$ESS =c(PS_MAP_ESS,MAP_ESS)
  
  
  if(getlimits==TRUE){
    limits=c(PS_MAP_fit[[3]]-0.25,PS_MAP_fit[[3]]+0.25)
    results=limits
  }
  
  return(results)
  
  
  
}





#two-arm trials
# results.trt = list()
# results.no = list()
# rejection.trt=NULL
# rejection.no=NULL
# tic()
# 
# 
# 
# for(i in 1:n_sims) {
#   results.trt[[i]] = RunSimulation(seed = i,trt.effect = 2,niter = niter)
#   #rejection.trt[i]=test_results( results.trt[[i]],cutoff=0.025)
# }
# 
# 
# for(i in 1:n_sims) {
#   results.no[[i]] = RunSimulation(seed = i,trt.effect = 0,niter = niter)
#   #rejection.no[i]=test_results( results.no[[i]],cutoff=0.025)
# }
# toc()
# 
# 
# 
# 
# 
# #mean(rejection.no)
# #mean(rejection.trt)
# rejection.trt = sapply(results.trt,FUN = test_results,cutoff = 0.025)
# rejection.no = sapply(results.no,FUN = test_results,cutoff = 0.025)
# 
# 
# data.frame(power = rowMeans(rejection.trt),
#            type1 = rowMeans(rejection.no), row.names = names(results.trt[[1]]))
# 
# 
# data.frame(power = apply(rejection.trt,MARGIN = 1,FUN = sum)/dim(rejection.trt)[2],
#            type1 = apply(rejection.no,MARGIN = 1,FUN = sum)/dim(rejection.no)[2], row.names = names(results.no[[1]]))
# 
