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

#Source functions for ps-power prior method
#source("ps-power_functions.R")
#source("ps-power_stan.R")


test_results = function(results, cutoff) {
  
  # Function to calculate if a posterior credible interval excludes zero
  # results: posterior as MCMC sample
  # cutoff: quantile value
  
  reject = rep(NA,length(results))
  for(k in 1:length(results)) {
    
    lower = summary(results[[k]][,'effect'], quantiles = cutoff)$quantiles
    upper = summary(results[[k]][,'effect'], quantiles = 1-cutoff)$quantiles
    reject[k] = !((lower < 0 ) & (upper > 0))
  }
  
  return(reject)
}


X_convert = function(X, binary_col) {
  
  # Convert continuous covariates to binary with cut point of 0
  # X: n*p matrix of covariates for n subjects, each has p covariates
  # binary_col: column indicator
  
  for (i in 1:length(binary_col)) {
    current_col = X[, i]
    X[which(current_col > 0), i] = 1
    X[which(current_col <= 0), i] = 0
  }
  return(t(X))
}




gen_Y = function(X, beta, var_eps,trt.effect=0,study.effect=0) {
  
  # Generate the binary outcome variable given the covariates
  # X: p*n matrix of covariates
  # beta: p*1 vector, coefficients for the covariates
  # var_eps: variance of the random noise
  # trt.effect: treatment effect
  # study.effect: study effect
  
  N = dim(X)[1]
  Y = matrix(0, N)
  for (i in 1:N) {
    xb = t(beta) %*% X[i, ]
    logit_P = xb +trt.effect+ study.effect+rnorm(1, 0, var_eps)
    P=inv.logit(logit_P)
    Y[i, ]=rbinom(1,1,P)
  }
 
  return(Y)
}



summarize_results = function(results, theta.true,cutoff,parameter) {
  
  # Summarize MCMC results and returns mean, bias, MSE, credible interval width, coverage
  # results: MCMC results as a list
  # theta.true: true effect
  # cutoff: cutoff quantile
  # parameter: name of parameter of interest 
  ## fixed to be 'effect' in a two-arm trial
  
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
    theta_MSE[k] = (mean(theta_draws[[1]])- theta.true) ^ 2
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




naive = function(response,n_patients, parameters = "effect", niter) {
  # Places a vague N(0,10000) prior on the logit of treatment and control response rate and returns MCMC object  
  
  # response: number of response from each trial as a vector, last two elements are control and trt from the current study
  # n_patients: number of patients for each trial as a vector,last two elements are # of patients in control and trt arm of the current study
  
  # parameters: vector of parameter names you want to monitor in the MCMC
  ## choices: 'effect': treatment effect on probability scale
  # 'beta0': log odds for control 
  # 'theta.trt' : log odds for treatment 
  # niter: number of MCMC iterations, will thin by 5 and burnin by 20% of niter
  
  n.trials = length(response)
  
  if(n.trials == 2) {
    dataTemp = list(
      "resp" = response,
      "n" = n_patients
    )
  } else {
    dataTemp = list(
      "resp" = c(sum(response[1:(n.trials - 1)]), response[n.trials]),
      "n" = c(sum(n_patients[1:(n.trials - 1)]),n_patients[n.trials])
    )
  }
  
  
  model <-
    jags.model(
      file = "./codes/PSMAP_Binary/PSMAP_Binary/naive.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  naive_model <-
    coda.samples(
      model,
      variable.names = parameters,
      thin = 5,
      n.iter = niter,
      progress.bar = 'none'
    )
  
  return(naive_model)
}



MAP.prior = function(response,n_patients, niter, std_heterogeneity = 1, parameters = 'p.pred') {
  
  # Obtain prior samples and ESS when using the MAP prior, returns the MCMC object contain the parameters
  # listed in "parameters"
  # response: number of response from each trial as a vector, last two elements are control and trt from the current study
  # n_patients: number of patients for each trial as a vector,last two elements are # of patients in control and trt arm of the current study
  # niter: number of MCMC iterations, will thin by 5 and burnin by 20% of niter
  # std_heterogeneity: hyper-parameter, std dev of half normal placed on tau
  # parameters: 'p.pred': a new sample from the prior 

  
    N=length(response)
  
  
  dataTemp = list(
    "resp" =response,
    "n" = n_patients,
    N = N
  )
  dataTemp$std_heterogeneity <- std_heterogeneity
  
  model <-
    jags.model(
      file = "./codes/PSMAP_Binary/PSMAP_Binary/MAP_Prior.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  MAP_model <-
    coda.samples(
      model,
      variable.names = parameters,
      thin = 5,
      n.iter = niter,
      progress.bar = 'none'
    )
  
  mix.res <- automixfit(MAP_model[[1]][,1], Nc=1:5, thresh=0, type="beta")
  ess.res=ess(mix.res, method="elir")
  return(list(MAP_model,ess.res))

}




MAP_fit = function(response,n_patients, niter,std_heterogeneity = 1, parameters = 'effect') {
  # Fits model using MAP prior on logit(p_control) and vague prior on logit(p_trt) and returns the MCMC object  
  # response: number of response from each trial as a vector, last two elements are control and trt from the current study
  # n_patients: number of patients for each trial as a vector,last two elements are # of patients in control and trt arm of the current study

  # niter: number of MCMC iterations, will thin by 5 and burnin by 20% of niter
  
  # std_heterogeneity: std for the half-Normal prior placed on tau
  
  # parameters: vector of parameters to return
  ## 'effect': treatment effect on probability scale
  # 'beta0': mean of the random effects model for controls
  # 'tau': precision of the random effects
  # 'delta': the random effects for each historical and current control 
  
  N=length(response)
  
  dataTemp = list(
    "resp" = response,
    "n" = n_patients,
    N = N
  )
  dataTemp$std_heterogeneity <- std_heterogeneity
  
  model <-
    jags.model(
      file = "./codes/PSMAP_Binary/PSMAP_Binary/MAP_Binary.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  MAP_model <-
    coda.samples(
      model,
      variable.names = parameters,
      thin = 5,
      n.iter = niter,
      progress.bar = 'none'
    )
  
  
  
  return(MAP_model)
}








power_fit = function(response, n_patients, niter, alpha = 0.1, a0 = 1, b0 = 1,parameters = 'effect') {
  
  # Fits model using power prior with vague normal prior on logit of trt and returns the MCMC object   
  # response: number of response from each trial as a vector, last two elements are control and trt from the current study
  # n_patients: number of patients for each trial as a vector,last two elements are # of patients in control and trt arm of the current study
  # niter: number of MCMC iterations, will thin by 5 and burnin by 20% of niter 
  # alpha: the power to raise the historic data to
  # a0/b0: the hyper-parameters for the beta distribution placed on original control response 
  # parameters: vector of parameter names you want to monitor in the MCMC
  ## choices: 
  # 'effect': treatment effect on probability scale
  # 'beta0': control log odds
  # 'theta.trt': treatment log odds
  # 'p[1]': control response rate
  # 'p[2]': treatment response rate
  

  n_hist_trials = length(response)-2
  
  dataTemp = list(
    "resp" = response[(n_hist_trials + 1):(n_hist_trials + 2)],
    "n" = n_patients[(n_hist_trials + 1):(n_hist_trials + 2)],
    "hist" = sum(response[1:n_hist_trials]),
    "n.hist" = sum(n_patients[1:n_hist_trials]),
    "alpha" = alpha,
    "a0" = a0,
    "b0" = b0
  )
  
  
  model <-
    jags.model(
      file = "./codes/PSMAP_Binary/PSMAP_Binary/PowerPrior.txt",
      data = dataTemp,
      n.chains = 1,
      n.adapt = 0.2 * niter,
      quiet = TRUE
    )
  update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
  power_model <-
    coda.samples(model,
                 variable.names = parameters,
                 thin = 5,
                 n.iter = niter, progress.bar = 'none')
  
  
  
  
  return(power_model)
  
}







PS_regroup<-function(X.hist.combine,Y.hist.combine,X.cur.con,Y.cur.con,p,S=5){
  
  # Estimate propensity scores and perform stratification, returns stratified data and stratum summary statistics
  # X.hist.combine: matrix of patient covariates, each row is a subject, historical studies
  # Y.hist.combine: vector of patient outcome,historical studies
  # X.cur.con: matrix of patient covariates, each row is a subject, current control arm
  # Y.cur.con: vector of patient outcome, current control arm
  # p: dimension of covariates
  
  #S=S
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
    X.hist.stratum [[i]] <- ana.ps$data[index,2:(p+1)]
    y.hist.stratum[[i]]<-Y.combine[index]
    
    indey2=which((ana.ps$data$group==1)&(ana.ps$data$"_strata_"==i))
    n.cur.stratum[i]=length(indey2)
    X.cur.stratum [[i]] <- ana.ps$data[indey2,2:(p+1)]
    y.cur.stratum[[i]]<-Y.combine[indey2]
  }
  
  stratum.ybar = sapply(y.hist.stratum,FUN = sum)
 
  stratum.ybar.cur = sapply(y.cur.stratum,FUN = sum)


  out=list(rS=rS,X.hist=X.hist.stratum, Y.hist=y.hist.stratum, Ybar.hist=stratum.ybar,
           Ybar.cur=stratum.ybar.cur, n.hist=n.hist.stratum, n.cur=n.cur.stratum,PS=ana.ps)
  return(out)
  
}







print("loading stan model...")

model <- stan_model(model_code = powerpsbinary, model_name = "powerpsbinary")

print("done")

PS_Power.fit=function(ps_data,outcomes,S=5,A){
  
  # Fits the PS-power prior with pre-specified sample size A and returns a posterior MCMC object 
  # ps_data: patient level covariates with propensity scores, all patients
  # outcomes: patient level outcome, all patients
  # S: number of strata
  # A: pre-specified number of patients intend to borrow
  
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
      type = "binary"
    )$post.theta
  
  
  
  results=theta_draws
  return(results)
}







PS_Power.fit2 = function(ybar.trt, n.trt, ps_data,outcomes,A, alpha0 = 1, beta0 = 1,
                         niter) {
  # Fits the PS-power prior and returns an MCMC object containing the posterior of treatment effect
  # ybar.trt/n.trt: the sample mean and sample size for treatement
  # ps_data: patient level covariates with propensity scores, all patients
  # outcomes: patient level outcome, all patients
  # A: pre-specified number of patents intend to borrow
  # alpha0/beta0: parameters for the beta prior on treatment mean
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  
  burn.in = (0.2)*niter
  # Create vectors to store parameters
  mu.trt = rbeta(niter+burn.in, ybar.trt+alpha0, n.trt-ybar.trt+beta0)
  
  
  mu.con=PS_Power.fit(ps_data,outcomes,S=5,A)
  mu.con= mu.con[seq(1,length(mu.con),by=5)]
  mu.trt=mu.trt[-(1:burn.in)]
  mu.trt= mu.trt[seq(1,length(mu.trt),by=5)]
  
  diff.effect = mu.trt -mu.con[1:length(mu.trt)]
  
  #Power.model = mcmc(data.frame(mu = mu.con[1:length(mu.trt)],mu.trt = mu.trt,effect = diff.effect))
  
  Power.model = mcmc(data.frame(effect = diff.effect))
  return(mcmc.list(Power.model))
}




PS_MAP.fit = function(tau.init, target.ESS, n.cur,ybar.hist, n.hist,overlap_coefficients, niter,lim, data.indiv) {
  
  # Obtain prior samples when using the PS-MAP prior with the specified target ESS
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # target.ESS: pre-specified target effective sample size

  # n.cur: number of subjects in each stratum from the current trial
  # ybar.hist: mean in each PS stratum from the historical trials
  # n.hist: number of subjects in each stratum from the historical trials
  # overlap_coefficients: overlap coefficients in each PS stratum
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # lim: range for value of tau.scale, as c(lower bound, upper bound)
  # data.indiv: outcome data for the current study
  
 
  n.strata = length(ybar.hist)
  total.cur=sum(n.cur)
  WS1=n.cur/total.cur
  ess.res=c(0)
  tau.res=c(0)
  direction=0

  tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.init #hyperparameter for HN(.)
  stop=FALSE
  low=lim[1]
  high=lim[2]
  tau.scale=tau.init
  while(stop==FALSE){
    theta.pred=list()
    # ess.res.str=rep(0,n.strata)
    
    #stop the loop if the ESS does not converge but the lower and upper bound are very close
    #which means too much variation in approximation using conjugate priors with small change in the hyper-parameter
    if((high-low<=0.01)){
      stop=TRUE
      tau.scale=tau.res[which.min(abs(ess.res-target.ESS))]
      tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.scale
    }
    
    
    #estimate stratum-specific MAP prior
    for(i in 1:n.strata){
      
      dataTemp = list(
        "resp" =ybar.hist[i],
        "n" = n.hist[i],
        N = 1
      )
      dataTemp$std_heterogeneity <- tau.prior[i]
      
      
      model <-
        jags.model(
          file = "./codes/PSMAP_Binary/PSMAP_Binary/MAP_Prior.txt",
          data = dataTemp,
          n.chains = 1,
          n.adapt = 0.2 * niter,
          quiet = TRUE
        )
      update(model, n.iter = 0.2 * niter, progress.bar = 'none') # burn in
      MAP_model <-
        coda.samples(model,
                     variable.names = "p.pred",
                     thin = 5,
                     n.iter = niter, progress.bar = 'none')
      
      theta.pred[[i]]=c(MAP_model[[1]][,1])
      
      #mix.res <- automixfit(theta.pred[[i]], Nc=1:4, thresh=0, type="beta")
      #ess.res.str[i]=ess( mix.res , method="elir")
      
    }
    
    
    theta.pred = do.call(cbind, theta.pred)
    theta=WS1%*%t(theta.pred)#overall prior as a weighted average of stratum-specific prior
    
    mix.res <- automixfit(theta[1,], Nc=1:5, thresh=0, type="beta")#approximation as mixture of beta
    
    ess.res=c(ess.res,ess(mix.res , method="elir"))#calculate ESS
    tau.res=c(tau.res,tau.scale)
    
    #Binary search for the value of tau.scale
    if(abs(ess.res[length(ess.res)]-target.ESS)<5){
      direction=0
      stop=TRUE
    }else if(ess.res[length(ess.res)]>target.ESS){
      direction=-1
      low=tau.scale
      tau.scale=(low+high)/2
      tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.scale
    }else{
      direction=1
      high=tau.scale
      tau.scale=(low+high)/2
      tau.prior=(min(overlap_coefficients)/overlap_coefficients)*tau.scale
    }
    
  }
  
  
  #calculate the posterior given prior and current data
  posterior.indiv<-postmix(mix.res, n=length(data.indiv),r=sum(data.indiv))
  results=list(mix.res,ess.res[length(ess.res)],posterior.indiv,tau.res[length(tau.res)])
  return(results)
}






PS_MAP.fit2=function(ybar.trt, n.trt, tau.init, target.ESS,  n.cur,ybar.hist,n.hist, overlap_coefficients,
                     lim, data.indiv, alpha0 = 1, beta0 = 1,niter) {
  # Fits the PS-MAP prior with target ESS and returns a posterior MCMC object of the treatment effect
  # ybar.trt/n.trt: the sample mean and sample size for treatment
  # tau.init: initial value for hyper-parameter of the half-Normal prior on tau^2
  # target.ESS: pre-specified target effective sample size
  # n.cur: number of subjects in each stratum from the current trial
  # ybar.hist: mean in each PS stratum from the historical trials
  # n.hist: number of subjects in each stratum from the historical trials
  # overlap_coefficients: overlap coefficients in each PS stratum
  # lim: range for value of tau.scale, as c(lower bound, upper bound)
  # data.indiv: outcome data for the current study
  # alpha0/beta0: parameters for the beta prior on treatment mean
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  
  
  burn.in = (0.2)*niter
  # Create vectors to store parameters
  mu.trt = rbeta(niter+burn.in, ybar.trt+alpha0, n.trt-ybar.trt+beta0)
  
  
  PS_MAP_fit= PS_MAP.fit(tau.init,target.ESS, n.cur, ybar.hist, n.hist, overlap_coefficients , 
                         niter,lim,data.indiv)
  
  PS_MAP_ESS=PS_MAP_fit[[2]] 
  PS_MAP_post=PS_MAP_fit[[3]]
  tau_scale_limits=PS_MAP_fit[[4]]
  mu.con=rmix(PS_MAP_post,n=niter)
  mu.con= mu.con[seq(1,length(mu.con),by=5)]
  mu.trt=mu.trt[-(1:burn.in)]
  mu.trt= mu.trt[seq(1,length(mu.trt),by=5)]
  diff.effect = mu.trt - mu.con[1:length(mu.trt)]
  
  PS_MAP.model = mcmc(data.frame(effect = diff.effect))
  
  fit.result=list(mcmc.list(PS_MAP.model),PS_MAP_ESS,tau_scale_limits)
  return(fit.result)
}




Two_arm_simulation = function(trt.effect=0, niter, MAP_ESS=90, Power_A=90, prior.std=0.5, alpha=0.1,limits=c(0.001,2),getlimits=FALSE) {
  
  # Run simulation in a two-arm setting and returns posterior of treatment effect based on different approaches
  # trt.effect: true value of treatment effect
  # niter: # of iterations burn in for 20% of niter and then thins by 5
  # MAP_ESS: target ESS for PS-MAP prior method
  # Power_A: target number of patients to borrow for PS-power prior method
  # prior.std: hyper parameter for prior on between-trial heterogeneity tau for MAP prior method
  # alpha: power parameter for Power prior method
  # limits: upper and lower values for tau.scale in PS-MAP prior method
  # getlimits: indicator for whether or not return the preliminary range for tau.scale in PS-MAP prior method
  
  ### Generate Data ----
  # Historic
  n.hist = 3 # number of historical trials
  n.hist.total=600 # total number of patients 
  n.hist.samples = rep(n.hist.total/n.hist,n.hist) # number of patients in each trial
  
  
  # Current 
  n.trt = 300 #number of patients in treatment
  n.con = 200 #number of patients in control
  
  #Generate covariates X
  p = 10 #dimension of the covariates
  rho = 0.1 #correlations of the covariates
  eps.var=1
  beta_0 = rep(0.1,p)
  
  mu.hist = cbind(rep(1,p),rep(1.4,p))
  S0 = diag(x = rep(1^2, p))
  S0[which(S0 == 0)] = rho * (1^2)
  
  
  mu.cur = matrix(1, p) # mean of covariates from the current study
  S1 = diag(x = rep(1^2, p))
  S1[which(S1 == 0)] = rho * (1^2)
  
  
  hist.list = list()
  X.hist =list()
  delta=c(0,0,0)
  
   for(i in 1:(n.hist-1)) {
     X.hist.init = mvrnorm(n.hist.samples[i], mu =mu.hist[,i], Sigma = S0)
     X.hist[[i]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
     hist.list[[i]] = gen_Y(X.hist[[i]], beta=beta_0, var_eps=eps.var,study.effect=delta[i])
   }
   
    mix.index=sum(rbinom(n.hist.samples[n.hist],1,1/3))
    X.hist.init_1 = mvrnorm(mix.index, mu =mu.hist[,n.hist-1], Sigma = S0)
    X.hist.init_2 = mvrnorm(n.hist.samples[n.hist]-mix.index, mu =mu.hist[,n.hist-2], Sigma = S0)
    X.hist.init=rbind(X.hist.init_1,X.hist.init_2)
    X.hist[[n.hist]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
    hist.list[[n.hist]] = gen_Y(X.hist[[n.hist]], beta=beta_0, var_eps=eps.var,study.effect=delta[n.hist])
  
  
  # for(i in 1:n.hist) {
  #   X.hist.init = mvrnorm(n.hist.samples[i], mu =mu.hist[,i], Sigma = S0)
  #   X.hist[[i]] = t(X_convert(X.hist.init, c(1, 2, 3, 4)))
  #   hist.list[[i]] = gen_Y(X.hist[[i]], beta=beta_0, var_eps=eps.var,study.effect=delta[i])
  # }
  # 
  
  X.hist.combine=abind(X.hist, along=1)
  Y.hist.combine=abind(hist.list, along=1)
  hist.ybar = sapply(hist.list,FUN = sum)

  
  
  # Current 
  #generate covariates for current study
  X.cur.init = mvrnorm((n.con+n.trt), mu =mu.cur, Sigma = S1)#assume trt and control from the same population
  X.cur = t(X_convert(X.cur.init, c(1, 2, 3, 4)))
  X.cur.con=X.cur[1:n.con,] #covaraites of the control group
  X.cur.trt=X.cur[(n.con+1):(n.trt+n.con),] #covariates of the treatment group
  
  
  #generate outcomes for control
  y.control = gen_Y(X.cur.con, beta=beta_0, var_eps=eps.var)
  ybar.control = sum(y.control)

  
  #generate outcomes for treatment
  y.treatment = gen_Y(X.cur.trt, beta=beta_0, var_eps=eps.var,trt.effect=trt.effect)
  ybar.trt = sum(y.treatment)

  
  
  #Calculate propensity scores
  regroup_data<-PS_regroup(X.hist.combine,Y.hist.combine,X.cur.con,y.control,p=p)
  
  
  
  
  ### Fit models ----
  
  pooled_model = naive(response=c(hist.ybar,ybar.control,ybar.trt),n_patients=c(n.hist.samples,n.con,n.trt), niter=niter,parameters = 'effect')
  
  
  current_model = naive(response=c(ybar.control,ybar.trt),n_patients=c(n.con,n.trt), niter=niter,parameters = 'effect')
  
 
  PS_MAP_fit=PS_MAP.fit2(ybar.trt = ybar.trt, n.trt=n.trt,
                         tau.init=1,target.ESS=MAP_ESS,
                         n.cur=regroup_data$n.cur, ybar.hist = regroup_data$Ybar.hist,n.hist=regroup_data$n.hist, 
                         overlap_coefficients = regroup_data$rS, 
                         lim=limits,data.indiv = y.control,
                         alpha0 = 1, beta0 = 1, niter = niter)
  
  PS_MAP_model=PS_MAP_fit[[1]]
  PS_MAP_ESS= PS_MAP_fit[[2]]
  
  MAP_model = MAP_fit(response=c(hist.ybar,ybar.control,ybar.trt),n_patients=c(n.hist.samples,n.con,n.trt), niter=niter,std_heterogeneity = prior.std, parameters = 'effect')
  
  
  MAP_ESS = MAP.prior(response= hist.ybar, n_patients=n.hist.samples,niter=niter, std_heterogeneity= prior.std, parameters = "p.pred")[[2]]
  
  
  
  PS_Power_model = PS_Power.fit2(ybar.trt=ybar.trt, n.trt=n.trt, ps_data=regroup_data$PS,outcomes=c(Y.hist.combine,y.control),A=Power_A, alpha0 = 1, beta0 = 1,
                           niter=niter)
  
  Power_model = power_fit(response=c(hist.ybar,ybar.control,ybar.trt),n_patients=c(n.hist.samples,n.con,n.trt),
                          niter=niter, alpha = alpha, a0 = 1, b0 = 1,parameters = 'effect')
  
  # Gather Results 
  
  results=list(Pooled=NULL, Current=NULL,PSMAP=NULL,MAP=NULL,PSPower=NULL,PowerPrior=NULL,ESS=NULL)
  
  results$Pooled=pooled_model
  results$Current=current_model
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


