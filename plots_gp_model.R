
## Check chain convergence
chains_diag <- lazymcmc::load_mcmc_chains(chainwd_gp, parTab,TRUE,1,mcmcPars_ct["adaptive_period"],
                                          multi=FALSE,chainNo=FALSE,PTchain = FALSE)
gelman_diagnostics(chains_diag$list)

## Load in main MCMC chains
chains <- lazymcmc::load_mcmc_chains(chainwd_gp, parTab,FALSE,1,mcmcPars_ct["adaptive_period"],
                                     multi=FALSE,chainNo=TRUE,PTchain = FALSE)
chain <- as.data.frame(chains$chain)
chain$sampno <- 1:nrow(chain)
chain_comb <- chain[chain$chain == 1,]
chain_comb$sampno <- 1:nrow(chain_comb)

## Plot fits to data
model_func <- create_posterior_func(parTab,obs_dat,NULL,inc_func_use,"model")
p_dist_fits <- plot_distribution_fits(chain_comb, obs_dat, model_func,100)
p_dist_fits

########################################
## Plot outputs from Gaussian process fit
########################################
## Get smoothed growth rates
samps <- sample(unique(chain_comb$sampno),n_samp)
trajs <- matrix(0, nrow=n_samp,ncol=length(times))
vl_trajs <- matrix(0, nrow=n_samp, ncol=length(ages))
detect_trajs <- matrix(0, nrow=n_samp, ncol=length(ages))

for(ii in seq_along(samps)){
  tmp_pars <- get_index_pars(chain_comb, samps[ii])
  vl <- viral_load_func(tmp_pars,ages,FALSE)
  vl_trajs[ii,]  <- extraDistr::rgumbel(length(vl),vl,tmp_pars["obs_sd"])
  detect_trajs[ii,] <- prop_detectable(ages, tmp_pars, vl)
  trajs[ii,] <- pmax(smooth.spline(inc_func_use(tmp_pars,times))$y,0.0000001)
}

trajs1 <- t(apply(trajs, 1, function(x) log(x[2:length(x)]/x[1:(length(x)-1)])))
trajs1_quants <- t(apply(trajs1, 2, function(x) quantile(x,c(0.025,0.5,0.975))))
trajs1_quants <- as.data.frame(trajs1_quants)
trajs1_quants$t <- 1:nrow(trajs1_quants)
colnames(trajs1_quants) <- c("lower","median","upper","t")
## Growth rate plot
p_gr <- ggplot(trajs1_quants) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) + 
  geom_line(aes(x=t,y=median),col=AAAS_palette["blue1"]) + 
  export_theme +
  geom_hline(yintercept=0,linetype="dashed") +
  ylab("Daily growth rate") +
  xlab("Time") +
  coord_cartesian(ylim=c(-0.5,0.5))

detect_trajs1 <- t(apply(detect_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
detect_trajs1 <- as.data.frame(detect_trajs1)
colnames(detect_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
detect_trajs1$t <- ages

p_detect <- ggplot(detect_trajs1) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  xlab("Days since infection") +
  ylab("Proportion detectable") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme

vl_trajs1 <- t(apply(vl_trajs, 2, function(x) quantile(x, c(0.025,0.25,0.5,0.75,0.975))))
vl_trajs1 <- as.data.frame(vl_trajs1)
colnames(vl_trajs1) <- c("lower","lower_mid","median","upper_mid","upper")
vl_trajs1$t <- ages

p_vl <- ggplot(vl_trajs1) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25) +
  geom_ribbon(aes(x=t,ymin=lower_mid,ymax=upper_mid),alpha=0.5) +
  geom_line(aes(x=t,y=median)) +
  xlab("Days since infection") +
  scale_y_continuous(trans="reverse") +
  scale_x_continuous(limits=c(0,50),breaks=seq(0,50,by=10)) +
  coord_cartesian(ylim=c(40,0)) +
  ylab("Ct value") +
  geom_vline(xintercept=5,linetype="dashed") +
  export_theme

trajs_quants <- t(apply(trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
trajs_quants <- as.data.frame(trajs_quants)
trajs_quants$t <- 1:nrow(trajs_quants)
colnames(trajs_quants) <- c("lower","mid_lower","median","mid_upper","upper","t")

obs_dat1 <- obs_dat %>% filter(ct < pars["intercept"])
obs_dat_all <- obs_dat1 %>% left_join(obs_dat1 %>% group_by(t) %>% tally())

p_dat <- ggplot(obs_dat_all) + 
  geom_violin(aes(x=t,group=t,y=ct,fill=n),scale="width",
              alpha=0.5,color="black",size=0.1,
              draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_dotplot(aes(x=t, y=ct,group=t),binaxis="y",
               binwidth=1,stackdir="center",binpositions="all",dotsize=0.1) +
  scale_y_continuous(trans="reverse",limits=c(42, 5),expand=c(0,0)) +
  scale_fill_gradient(low=AAAS_palette["blue1"],high=AAAS_palette["red1"]) +
  geom_hline(yintercept=40,linetype="dashed") +
  export_theme +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = c(0.05,0.5)) +
  xlab("Time of sample") +
  ylab("Ct value") +
  labs(tag="B")
p_dat

trajs_melted <- reshape2::melt(trajs)
colnames(trajs_melted) <- c("samp","t","inc")
trajs_melted$t <- times[trajs_melted$t]
## Growth rate plot
p_inc <- ggplot(trajs_quants) + 
  geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.25,fill=AAAS_palette["blue1"]) + 
  geom_ribbon(aes(x=t,ymin=mid_lower,ymax=mid_upper),alpha=0.5,fill=AAAS_palette["blue1"]) + 
  geom_line(aes(x=t,y=median),col=AAAS_palette["blue1"]) + 
  export_theme +
  ylab("Relative probability of infection") +
  xlab("Time") +
  scale_y_continuous(expand=c(0,0))  +
  scale_x_continuous(limits=range(seir_dynamics$seir_outputs$step)) +
  labs(tag="C")
