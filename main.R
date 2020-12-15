########################################
## SCRIPT TO SIMULATE LINE LIST DATA ON CLUSTER
## 11th December 2020
## ------ Aims ---------
#' 1. Loads the virosolver R package
#' 2. Read in the correct parTab data frame, which controls the incidence and viral load simulations
#' 3. Simulate a full line list data set, giving infection times etc for each individual in the population
#' 4. Simulate Ct values for the sub-setted line list, giving a Ct value for each observed individual
########################################
########################################
## 1. Headers
########################################
library(tidyverse)
library(ggplot2)
library(extraDistr)
library(lazymcmc) ## install from devtools::install_github("jameshay218/lazymcmc")
library(patchwork)
library(ggthemes)
library(odin) ## install from CRAN
library(doParallel)
library(fitdistrplus)

HOME_WD <- "~/Documents/GitHub/"
devtools::load_all(paste0(HOME_WD,"/virosolver"))

## Where to perform the simulations
main_wd <- paste0(HOME_WD,"/virosolver_ons")
setwd(main_wd)

## What to call this run
task_id <- 1
run_name <- paste0("simtest_",task_id)

## Where to save the simulated data
save_wd <- paste0(HOME_WD, "/virosolver_paper/data/ONS_sim/")
save_file <- paste0(save_wd, run_name)
if(!file.exists(save_wd)) dir.create(save_wd,recursive = TRUE)

## Where to save MCMC chains and outputs
chainwd_gp <- paste0(main_wd,"/mcmc_chains/",run_name,"_gp")
chainwd_exp <- paste0(main_wd,"/mcmc_chains/",run_name,"_exp")
plot_wd <- paste0(main_wd,"/plots/",run_name)
if(!file.exists(chainwd_gp)) dir.create(chainwd_gp,recursive = TRUE)
if(!file.exists(chainwd_exp)) dir.create(chainwd_exp,recursive = TRUE)
if(!file.exists(plot_wd)) dir.create(plot_wd,recursive = TRUE)

## Arguments for this run
set.seed(1)

##################################################
## NOTE THIS - change flag to FALSE if want to include undetectable Cts
use_pos <- FALSE
n_samp <- 1000 ## How many posterior samples to use for plots etc
rerun_exp <- FALSE
rerun_gp <- TRUE

##################################################
## NOTE THIS BIT TO USE FOREACH LOOP
## Manage MCMC runs and parallel runs
n_clusters <- 10
cl <- parallel::makeCluster(n_clusters, setup_strategy = "sequential")
registerDoParallel(cl)
##################################################


## Load functions for line list simulation
source("code/linelist_sim_funcs.R")
source("code/odin_funcs.R")
source("code/plot_funcs.R")
## Priors for all models - EDIT THIS FILE TO CHANGE PRIORS!
source("code/priors.R")

########################################
## 2. Model parameters and simulation settings
########################################
## Model parameters
model_pars <- read.csv("pars/partab_seir_model.csv")
pars <- model_pars$values
names(pars) <- model_pars$names

## Simulation parameters
population_n <- 1000000
## Over the course of 200 days
times <- 0:250

## Sampling procedure
sampling_frequency <- 14
sampling_number <- 2000

########################################
## IMPORTANT CHECK
## - Check that the assumed viral kinetics are in line
##   with your data. This means things like peak Ct value,
##   waning rate, level of the "plateau" phase, and the 
##   limit of detection
########################################
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

########################################
## 3. Full simulated line list
########################################
source("simulate_data.R")
head(obs_dat)
obs_dat <- obs_dat %>% filter(t < 150)
ages <- 1:max(obs_dat$t)
times <- 0:max(obs_dat$t)

########################################
## 4. Fit single cross sections to simulated data
########################################
#source("fit_exp_model.R")
#source("plots_exp_model.R")
#p_exp_betas
#ggsave("plots/exp_betas_example.png",p_exp_betas,width=8,height=5,units="in",dpi=300)

########################################
## 5. Fit full Gaussian process model to simulated data
########################################
source("fit_gp_model.R")
source("plots_gp_model.R")
p_dat1 <- p_dat + scale_x_continuous(limits=c(0,250))
p_vl
p_detect
p_gr
seir_dynamics$plot/p_inc

(p_dat + scale_x_continuous(limits=c(0,250)))/p_inc

ggsave("plots/gp_example.png",seir_dynamics$plot/p_inc,width=8,height=10,units="in",dpi=300)
