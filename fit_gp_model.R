ages <- 1:max(obs_dat$t)
times <- 0:(max(obs_dat$t))


nchains <- 3
## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=10000,"popt"=0.44,"opt_freq"=2000,
                 "thin"=10,"adaptive_period"=10000,"save_block"=100)

## Model parameters and simulation settings
inc_func_use <- gaussian_process_model
prior_func_use <- prior_func_hinge_gp

## GP model parameters for fitting
parTab <- read.csv(paste0(main_wd,"/pars/partab_gp_model.csv"))
if(!use_pos) {
  parTab[parTab$names == "overall_prob","fixed"] <- 0
}
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

## This is for the GP version
mat <- matrix(rep(times, each=length(times)),ncol=length(times))
t_dist <- abs(apply(mat, 2, function(x) x-times))
t_dist[upper.tri(t_dist)] <- Inf
parTab <- bind_rows(parTab[parTab$names != "prob",], parTab[parTab$names == "prob",][1:length(times),])
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

inc_func_use <- gaussian_process_model
## Check that posterior function solves correctly
f <- create_posterior_func(parTab, obs_dat, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=use_pos,
                           t_dist=t_dist)
f(parTab$values)

pars1 <- pars
true_probs <- seir_dynamics$incidence/population_n
true_probs <- true_probs[which(names(pars1) == "prob")]
pars1[which(names(pars1) == "prob")] <- reverse_gp_model(true_probs, pars1, times)
f(pars1)

## Run for each chain
chains <- NULL
if(rerun_gp){
  res <- foreach(j=1:nchains,.packages = c("lazymcmc","extraDistr","tidyverse","patchwork")) %dopar% {
    devtools::load_all(paste0(HOME_WD,"/virosolver"))
    
    ## Get random starting values
    startTab <- generate_viable_start_pars(parTab,obs_dat,
                                           create_posterior_func,
                                           inc_func_use,
                                           prior_func_use,
                                           t_dist=t_dist,
                                           use_pos=use_pos)
    covMat <- diag(nrow(startTab))
    mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
    
    output <- run_MCMC(parTab=startTab,
                       data=obs_dat,
                       INCIDENCE_FUNC=inc_func_use,
                       PRIOR_FUNC = prior_func_use,
                       solve_likelihood=TRUE,
                       mcmcPars=mcmcPars_ct,
                       filename=paste0(chainwd_gp,"/",run_name,"_chainno_",j),
                       CREATE_POSTERIOR_FUNC=create_posterior_func,
                       mvrPars=NULL,
                       OPT_TUNING=0.2,
                       use_pos=use_pos,
                       t_dist=t_dist)
    
    ## Read in chain and remove burn in period
    chain <- read.csv(output$file)
    chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
    chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
    chains[[j]] <- chain
    chain <- do.call("bind_rows",chains)
  } 
}
