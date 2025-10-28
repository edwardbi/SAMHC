
source("PSMAP_Normal_functions.R")


#Specify simulation settings & design parameters
true.trt.effect=0
ESS=100
prior.std=0.1
alpha=0.167
n_sims = 10
seed = 110





set.seed(seed)
niter = 100000
n.methods=6
result_final=list()
final_means=matrix(0,n_sims,n.methods)
colnames(final_means) <- c('Pooled', 'NoBorrow','PS-MAP', 'MAP', 'PS-Power', 'Power')
final_bias=matrix(0,n_sims,n.methods)
colnames(final_bias) <- c('Pooled', 'NoBorrow','PS-MAP', 'MAP', 'PS-Power', 'Power')
final_MSE=matrix(0,n_sims,n.methods)
colnames(final_MSE) <- c('Pooled', 'NoBorrow','PS-MAP', 'MAP', 'PS-Power', 'Power')
final_width=matrix(0,n_sims,n.methods)
colnames(final_width) <- c('Pooled', 'NoBorrow','PS-MAP', 'MAP', 'PS-Power', 'Power')
final_coverage=matrix(0,n_sims,n.methods)
colnames(final_coverage) <- c('Pooled', 'NoBorrow','PS-MAP', 'MAP', 'PS-Power', 'Power')
final_rejection=matrix(0,n_sims,n.methods)
colnames(final_rejection) <- c('Pooled', 'NoBorrow','PS-MAP', 'MAP', 'PS-Power', 'Power')
final_ESS=matrix(0,n_sims,n.methods)
colnames(final_ESS) <- c('PS-MAP', 'MAP','PS-MAP', 'MAP','PS-MAP', 'MAP')



tic()
#get the approximated range for tau.scale in PS-MAP method
limits = Two_arm_simulation(trt.effect=true.trt.effect, niter=niter, MAP_ESS=ESS, Power_A=ESS, prior.std=prior.std, alpha=alpha,limits=c(0.00001,2),getlimits=TRUE) 

#Run simulations
for(i in 1:n_sims) {
  sumresults=matrix(0,nrow=6,ncol=n.methods)
  results_i = Two_arm_simulation(trt.effect=true.trt.effect, niter=niter, MAP_ESS=ESS, Power_A=ESS, prior.std=prior.std, alpha=alpha,limits=limits) 
  result_final[[i]]=results_i
  sumresults_i=summarize_results(results_i,theta.true = true.trt.effect,cutoff=0.025,parameter="effect")
  final_means[i,]=as.numeric(sumresults_i$mean)
  final_bias[i,]=as.numeric(sumresults_i$bias)
  final_MSE[i,]=as.numeric(sumresults_i$MSE)
  final_width[i,]=as.numeric(sumresults_i$CI_width)
  final_coverage[i,]=as.numeric(sumresults_i$CI_coverage)
  final_rejection[i,]=as.numeric(sumresults_i$reject)
  final_ESS[i,]=results_i$ESS
  print(i)
}

toc()



#See Results
colMeans(final_means)
colMeans(final_bias)*100
colMeans(final_MSE)*100
colMeans(final_rejection)*100
colMeans(final_ESS)

