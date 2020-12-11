nchains <- 3
## MCMC parameters for Ct model fits
mcmcPars_ct_exp <- c("iterations"=50000,"popt"=0.44,"opt_freq"=2000,
                     "thin"=10,"adaptive_period"=20000,"save_block"=100)
max_age <- 35
inc_func_use <- exponential_growth_model
prior_func_use <- prior_func_hinge_exp

## GP model parameters for fitting
parTab <- read.csv("pars/partab_exp_model.csv")
if(!use_pos) {
  parTab[parTab$names == "overall_prob","fixed"] <- 0
}
pars <- parTab$values
names(pars) <- parTab$names

## Means for priors
means <- parTab$values
names(means) <- parTab$names

obs_times <- unique(obs_dat$t)
ages <- 1:max(obs_dat$t)
times <- 0:max(obs_dat$t)

## Check that posterior function solves correctly
f <- create_posterior_func(parTab, obs_dat, prior_func_use, 
                           inc_func_use,solve_ver="likelihood",
                           use_pos=use_pos,
                           t_dist=t_dist)
f(pars)
if(rerun_exp){
  res <- foreach(i=seq_along(obs_times),.packages = c("extraDistr","tidyverse","patchwork","lazymcmc")) %dopar% {
    ## Need to read in R package
    devtools::load_all(paste0(HOME_WD,"/virosolver"))
    timepoint <- obs_times[i]
    runname_use <- paste0(run_name,"_time_",timepoint)
    dir.create(paste0(chainwd_exp,"/",timepoint),recursive = TRUE)
    obs_dat_use <- obs_dat %>% filter(t == timepoint)
    
    ## Observation times
    if(!is.na(max_age)){
      obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t), t = t + max_age)
    }
    ages <- 1:max(obs_dat_use$t)
    times <- 0:max(obs_dat_use$t)
    parTab[parTab$names == "t0",c("upper_bound","upper_start")] <- min(obs_dat_use$t)
    
    chains <- NULL
    for(j in 1:nchains){
      ## Get random starting values
      startTab <- generate_viable_start_pars(parTab,obs_dat_use,
                                             create_posterior_func,
                                             inc_func_use,
                                             prior_func_use,
                                             t_dist=NULL,
                                             use_pos=use_pos)
      
      covMat <- diag(nrow(startTab))
      mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
      output <- run_MCMC(parTab=startTab,
                         data=obs_dat_use,
                         INCIDENCE_FUNC=inc_func_use,
                         PRIOR_FUNC = prior_func_use,
                         solve_likelihood=TRUE,
                         mcmcPars=mcmcPars_ct_exp,
                         filename=paste0(chainwd_exp,"/",timepoint,"/",runname_use,"_chainno_",j),
                         CREATE_POSTERIOR_FUNC=create_posterior_func,
                         mvrPars=mvrPars,
                         OPT_TUNING=0.2,
                         use_pos=use_pos,
                         t_dist=NULL)
      
      ## Read in chain and remove burn in period
      chain <- read.csv(output$file)
      chain <- chain[chain$sampno > mcmcPars_ct["adaptive_period"],]
      chain$sampno <-chain$sampno + max(chain$sampno)*(j-1)
      chain$chain <- j
      chains[[j]] <- chain
    }
    chain <- do.call("bind_rows",chains)
  }
} 

