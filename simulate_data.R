
## Simulate SEIR dynamics, incidence, growth rates, Rt etc
## Use "ode" for deterministic
## Use "odin" for stochastic
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=pars, ver="odin",switch_model=FALSE)

## Simulate onset times, confirmation delays etc
## This returns a tibble with line list entries for **every** individual in the population
complete_linelist <- virosolver::simulate_observations_wrapper(seir_dynamics$incidence,times=times,
                                                               population_n=population_n)

########################################
## Simulate observation process
########################################
## This bit of code generates linelist data for some specified observation process.
## Here, we simulate sampling 1000 people at random from the population every 14 days
## arguments changed above
sample_probs <- c(rep(0, sampling_frequency-1),sampling_number/population_n)
sample_probs <- rep(sample_probs, length(times)/sampling_frequency +1)
sample_probs <- sample_probs[1:length(times)]
frac_report <- tibble(t=times,prob=sample_probs)
frac_report <- frac_report %>% filter(t >= 50 & t <= 160)

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