#!/usr/bin/env Rscript

library(argparser)
library(magrittr)
source("PSMAP_Binary_functions.R")


p = arg_parser("Program description") %>%
  add_argument("true.trt.effect", help="trt effect in current trial", type="numeric") %>%
  add_argument("ESS", help="Target effective sample size", type="numeric") %>%
  add_argument("prior.std", help="hyper parameter in MAP Prior", type="numeric") %>%
  add_argument("alpha", help="scale parameter in Power prior", type="numeric") %>%
  add_argument("part", help="Part number", type="numeric") 



argv = parse_args(p)

task_id = Sys.getenv("SLURM_ARRAY_TASK_ID") %>% as.numeric
part=argv$part

seed = 10*task_id + part

set.seed(seed)
n_sims = 50
niter = 100000
n.methods=6
true.trt.effect=argv$true.trt.effect
ESS=argv$ESS
prior.std=argv$prior.std
alpha=argv$alpha

result_final=list()
final_means=matrix(0,n_sims,n.methods)
final_bias=matrix(0,n_sims,n.methods)
final_MSE=matrix(0,n_sims,n.methods)
final_width=matrix(0,n_sims,n.methods)
final_coverage=matrix(0,n_sims,n.methods)
final_rejection=matrix(0,n_sims,n.methods)
final_ESS=matrix(0,n_sims,n.methods)



tic()

limits = Two_arm_simulation(trt.effect=true.trt.effect, niter=niter, MAP_ESS=ESS, Power_A=ESS, prior.std=prior.std, alpha=alpha,limits=c(0.00001,2),getlimits=TRUE) 
print(limits)
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





result_sum=list(final_means, final_bias,final_MSE,final_width,final_coverage, final_rejection,final_ESS, result_final)
filename=paste0("results_",true.trt.effect,"_",ESS,"_",seed,".RData")
save(result_sum, file =filename)


