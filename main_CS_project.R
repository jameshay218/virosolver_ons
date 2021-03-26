########################################
## SCRIPT TO SIMULATE LINE LIST DATA AND CT VALUES 
## FROM AN SEIR PROCESS, THEN RE-ESTIMATE Rt
## 26th March 2021
## ------ Structure ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame, which controls the incidence curve and viral load simulations
#' 3. Simulate a full line list data set, giving infection times etc for each individual in the population
#' 4. Simulate Ct values for the sub-setted line list, giving a Ct value for each observed individual
#' 5. Feed these Ct values to the inference framework and re-estimate the simulated epidemic trajectory

########################################
## 1. Headers
########################################
## Install these packages
library(tidyverse)
library(ggplot2)
library(extraDistr)
library(patchwork)
library(ggthemes)
library(odin) ## install from CRAN
library(doParallel)
library(fitdistrplus)

library(lazymcmc) ## install from devtools::install_github("jameshay218/lazymcmc")
library(virosolver) ## install from devtools::install_github("jameshay218/virosolver")

HOME_WD <- "~/Documents/GitHub/"

## Where to perform the simulations
main_wd <- paste0(HOME_WD,"/virosolver_ons")
setwd(main_wd)

## NOTE: no need to open these files unless you really want to.
## Just trust that they provide functions for your analyses.
## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
## Probably fine if you just leave these unchanged.
source("code/priors.R")

## What to call this run
task_id <- 1
run_name <- paste0("CS_project_",task_id)

## Where to save the simulated data
save_wd <- paste0(HOME_WD, "/virosolver_paper/data/ONS_sim/")
save_file <- paste0(save_wd, run_name)
if(!file.exists(save_wd)) dir.create(save_wd,recursive = TRUE)

## Where to save MCMC chains and outputs
chainwd_seir <- paste0(main_wd,"/mcmc_chains/",run_name,"_seir")
chainwd_exp <- paste0(main_wd,"/mcmc_chains/",run_name,"_exp")
plot_wd <- paste0(main_wd,"/plots/",run_name)
if(!file.exists(chainwd_seir)) dir.create(chainwd_seir,recursive = TRUE)
if(!file.exists(chainwd_exp)) dir.create(chainwd_exp,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

## Arguments for this run
set.seed(1)

## NOTE THIS - change flag to FALSE if want to include undetectable Cts
use_pos <- FALSE
n_samp <- 1000 ## How many posterior samples to use for plots etc
rerun_seir <- TRUE

## NOTE THIS BIT TO USE FOREACH LOOP
## Manage MCMC runs and parallel runs
n_clusters <- 10
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
model_pars <- read.csv("pars/partab_seir_model.csv")
View(model_pars)
pars <- model_pars$values
names(pars) <- model_pars$names

## IMPORTANT CHECK
## - Check that the assumed viral kinetics are in line
##   with your data. This means things like peak Ct value,
##   waning rate, level of the "plateau" phase, and the 
##   limit of detection
test_ages <- seq(1,50,by=1)
cts <- viral_load_func(pars, test_ages)
prop_detect <- prop_detectable(test_ages,pars, cts)
p1 <- ggplot(data.frame(ct=cts,t=test_ages)) + 
  geom_line(aes(x=t,y=ct)) + 
  scale_y_continuous(trans="reverse",
                     limits=c(40,0)) +
  ylab("Modal Ct value") +
  xlab("Days since infection")
p2 <- ggplot(data.frame(p=prop_detect,t=test_ages)) + 
  geom_line(aes(x=t,y=p)) + 
  ylab("Proportion of infections still detectable") +
  xlab("Days since infection")
p1/p2

## Simulation parameters
population_n <- 1000000
## Over the course of 250 days
times <- 0:250

## Sampling procedure - how often do we take samples from the population?
sampling_frequency <- 14
## How many samples do we take on each sample day?
sampling_number <- 2000

########################################
## 3. Full simulated line list
########################################
## Simulate SEIR dynamics, incidence, growth rates, Rt etc
## Use "ode" for deterministic
## Use "odin" for stochastic model
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=pars, ver="odin",switch_model=FALSE)

## Look at simulated dynamics over time
seir_dynamics$plot

## Basic plot, showing the proportion of the population infected each day.
## This is the epidemic incidence curve that you're trying to estimate.
## So you can substitute this object for **any** incidence curve that
## you want to estimate.
plot(seir_dynamics$incidence/population_n)

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
## These entries give the **true** infection (or lack of infection) timings for each individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                               population_n=population_n)

########################################
## 4. Simulate observation process
########################################
## This bit of code generates linelist data for some specified observation process.
## Here, we simulate sampling 2000 people at random from the population every 14 days
## arguments changed above
sample_probs <- c(rep(0, sampling_frequency-1),sampling_number/population_n)
sample_probs <- rep(sample_probs, length(times)/sampling_frequency +1)
sample_probs <- sample_probs[1:length(times)]
frac_report <- tibble(t=times,prob=sample_probs)
frac_report <- frac_report %>% filter(t >= 50 & t <= 160)

## frac_report is a table, giving the proportion (prob) of the population
## sampled on day t
head(frac_report)

## This function takes the complete linelist and sub-samples a proportion
## of the entire population on the days specified by `frac_report`
observed_linelist <- simulate_reporting(complete_linelist, 
                                        frac_report=NULL,
                                        timevarying_prob=frac_report,
                                        solve_times=times, 
                                        symptomatic=FALSE)

## Simulate viral load/Ct value for each person in the line list
simulated_viral_loads <- simulate_viral_loads_wrapper(observed_linelist$sampled_individuals,
                                                      kinetics_pars=pars)

## Clean data to form expected by virosolver
obs_dat <- simulated_viral_loads %>% dplyr::select(sampled_time, ct_obs) %>%
  rename(t = sampled_time, ct=ct_obs) %>% arrange(t)

## Save simulated line list, Cts and SEIR dynamics
write_csv(seir_dynamics$seir_outputs, path=paste0(save_file,"_seir_outputs.csv"))
write_csv(complete_linelist, path=paste0(save_file,"_full_linelist.csv"))
write_csv(obs_dat, path=paste0(save_file,"_cts.csv"))

## Save SEIR plots
## Ct distribution plot
p_dat <- ggplot(obs_dat %>% filter(ct < pars["intercept"])) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(aes(x=t,y=ct),size=0.1,width=2,height=0) + 
  scale_y_continuous(trans="reverse") +
  export_theme +
  scale_x_continuous(limits=c(min(times),max(times)+50)) +
  ylab("Ct value") +
  xlab("Observation time")

ggsave(paste0(save_file,"_ct_dist.png"),p_dat,width=7,height=4,dpi=150)
ggsave(paste0(save_file,"_seir_model.png"),seir_dynamics$plot,width=7,height=7,dpi=150)
ggsave(paste0(save_file,"_growth_rates.png"),seir_dynamics$growth_rate_p,width=7,height=6,dpi=150)

head(obs_dat)
obs_dat <- obs_dat %>% filter(t < 150)
ages <- 1:max(obs_dat$t)
times <- 0:max(obs_dat$t)


########################################
## 5. Fit the SEIR model to the simulated data
########################################
## REMEMBER, obs_dat is your simulated dataset.
## Your aim is to use this dataset to re-estimate the
## contents of seir_dynamics.
nchains <- 3
## MCMC parameters for Ct model fits
mcmcPars_ct <- c("iterations"=50000,"popt"=0.44,"opt_freq"=2000,
                     "thin"=10,"adaptive_period"=20000,"save_block"=100)
inc_func_use <- solveSEIRModel_rlsoda_wrapper
prior_func_use <- prior_func_hinge_seir

## GP model parameters for fitting
parTab <- read.csv("pars/partab_seir_model.csv")
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
if(rerun_seir){
  res <- foreach(i=seq_along(obs_times),.packages = c("extraDistr","tidyverse","patchwork","lazymcmc","virosolver")) %dopar% {
    ## Need to read in R package
    timepoint <- obs_times[i]
    runname_use <- paste0(run_name,"_time_",timepoint)
    dir.create(paste0(chainwd_seir,"/",timepoint),recursive = TRUE)
    obs_dat_use <- obs_dat %>% filter(t == timepoint)
    
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
                         mcmcPars=mcmcPars_ct,
                         filename=paste0(chainwd_seir,"/",timepoint,"/",runname_use,"_chainno_",j),
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

########################################
## 6. Generate plots from the estimated posterior
########################################
res <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  chainwd_tmp <- paste0(chainwd_seir,"/",timepoint)
  chain <- lazymcmc::load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                      multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain <- as.data.frame(chain)
  res[[i]] <- chain
}

## Go through each time point that we estimated dynamics from
## Get the MCMC chains ie. the posteriors and generate some plots
## and analyses
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use <- runname_use <- paste0(run_name,"_time_",timepoint)
  
  obs_dat_tmp <- obs_dat_use <- obs_dat %>% filter(t == timepoint)
  
  ages <- 1:max(obs_dat_use$t)
  times <- 0:max(obs_dat_use$t)
  
  chain <- res[[i]]
  chain_comb <- chain
  chain_comb$sampno <- 1:nrow(chain_comb)
  chain1 <- chain
  chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
  
  
  ## Check convergence of MCMC chains. 
  ## We want to see lots of hairy caterpillars.
  ## Some of the parameters will be flat lines - these
  ## were the parameters that we fixed.
  p_trace <- chain1[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
    mutate(chain = as.factor(chain)) %>%
    pivot_longer(-c(sampno,chain)) %>%
    ggplot() +
    geom_line(aes(x=sampno,y=value,col=chain)) +
    facet_wrap(~name,scales="free_y")+
    scale_x_continuous(breaks=seq(min(chain$sampno),max(chain$sampno),length.out=5)) +
    export_theme
  
  ## Get smoothed growth rates
  samps <- sample(unique(chain_comb$sampno),n_samp)
  trajs <- matrix(0, nrow=n_samp,ncol=length(times))
  for(ii in seq_along(samps)){
    trajs[ii,] <- pmax(inc_func_use(get_index_pars(chain_comb, samps[ii]),times),0.0000001)
  }
  
  trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
  trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  trajs1_quants <- as.data.frame(trajs1_quants)
  trajs1_quants$t <- 1:nrow(trajs1_quants)
  colnames(trajs1_quants) <- c("lower","median","upper","t")
  
  ## Growth rate plot
  p_gr <- ggplot(trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median)) + 
    coord_cartesian(ylim=c(-0.5,0.5))
  
  trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
  trajs_quants <- as.data.frame(trajs_quants)
  trajs_quants$t <- 1:nrow(trajs_quants)
  colnames(trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")
  
  ## Incidence plot
  p_inc <- ggplot(trajs_quants) + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_ribbon(aes(x=t,ymin=mid_lower,ymax=mid_upper),alpha=0.5) + 
    geom_line(aes(x=t,y=median)) + 
    geom_line(data=tibble(x=times,y=seir_dynamics$incidence[seq_along(times)]/population_n),aes(x=x,y=y)) +
    export_theme +
    ylab("Per capita incidence") +
    xlab("Days since start") +
    coord_cartesian(ylim=c(0,0.03))
  
  vl_trajs <-  matrix(0, nrow=n_samp,ncol=length(ages))
  for(ii in 1:n_samp){
    tmp_pars <- get_index_pars(chain_comb, samps[ii])
    tmp <- viral_load_func(tmp_pars,ages,FALSE)
    tmp1 <- extraDistr::rgumbel(length(tmp),tmp, tmp_pars["obs_sd"])
    vl_trajs[ii,] <- tmp1
  }
  vl_trajs1_quants <- t(apply(vl_trajs, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
  vl_trajs1_quants <- as.data.frame(vl_trajs1_quants)
  vl_trajs1_quants$t <- 1:nrow(vl_trajs1_quants)
  colnames(vl_trajs1_quants) <- c("lower","median","upper","t")
  
  ## Growth rate plot
  p_vl <- ggplot(vl_trajs1_quants) + geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_line(aes(x=t,y=median))
  
  dir.create(paste0(plot_wd,"/traces/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/predictions/"),recursive = TRUE)
  dir.create(paste0(plot_wd,"/grs/"),recursive = TRUE)
  
  ggsave(paste0(plot_wd,"/traces/",runname_use,"_trace.png"),p_trace,width=7,height=4)
  ggsave(paste0(plot_wd,"/predictions/",runname_use,"_predictions.png"),p_dat/p_inc,width=7,height=7)
  ggsave(paste0(plot_wd,"/grs/",runname_use,"_grs.png"),p_gr,width=7,height=4)
}
