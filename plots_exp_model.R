
res <- NULL
beta_ests <- NULL
for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  chainwd_tmp <- paste0(chainwd_exp,"/",timepoint)
  chain <- lazymcmc::load_mcmc_chains(chainwd_tmp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                      multi=TRUE,chainNo=TRUE,PTchain = TRUE)$chain
  chain <- as.data.frame(chain)
  beta_ests[[i]] <- tibble(beta=chain$beta,t=timepoint)
  res[[i]] <- chain
}
beta_ests_combined <- do.call("bind_rows",beta_ests)

for(i in seq_along(obs_times)){
  timepoint <- obs_times[i]
  runname_use <- runname_use <- paste0(run_name,"_time_",timepoint)
  
  obs_dat_tmp <- obs_dat_use <- obs_dat1 %>% filter(t == timepoint)
  
  ## Observation times
  if(!is.na(max_age)){
    obs_dat_use <- obs_dat_use %>% mutate(t = t - min(t), t = t + max_age)
  }
  
  ages <- 1:max(obs_dat_use$t)
  times <- 0:max(obs_dat_use$t)
  
  chain <- res[[i]]
  chain_comb <- chain
  chain_comb$sampno <- 1:nrow(chain_comb)
  chain1 <- chain
  chain_comb <- chain_comb[,colnames(chain_comb) != "chain"]
  
  
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
  
  ## Growth rate plot
  p_inc <- ggplot(trajs_quants) + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) + 
    geom_ribbon(aes(x=t,ymin=mid_lower,ymax=mid_upper),alpha=0.5) + 
    geom_line(aes(x=t,y=median)) + 
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

p_exp_betas <- beta_ests_combined %>% 
  ggplot() + geom_violin(aes(x=t,y=beta,group=t),draw_quantiles=c(0.025,0.5,0.975),fill="grey70") +
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("35-day growth rate using exponential model") +
  xlab("Time") +
  export_theme


