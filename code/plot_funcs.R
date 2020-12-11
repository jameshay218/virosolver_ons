
## Create an epiweek calendar
dates <- seq(as.Date("2020-01-01"),as.Date("2020-12-31"),by="1 day")
epiweeks <- lubridate::epiweek(dates)
epi_calendar <- tibble(date=dates,week=epiweeks)
epi_calendar <- epi_calendar %>% group_by(week) %>% mutate(first_day=min(date))


export_theme <- theme_tufte() +
  theme(
    axis.text.x = element_text(size=7,family="sans"),
    axis.text.y=element_text(size=7,family="sans"),
    axis.title.x=element_text(size=8,family="sans",vjust=-1),
    axis.title.y=element_text(size=8,family="sans"),

    ## Axis lines
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),

    ## Title
    plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
    plot.tag = element_text(family="sans",size=10),

    ## Legends
    legend.title=element_text(size=8,family="sans",face="italic"),
    legend.text=element_text(size=8,family="sans"),
    legend.key.size= unit(0.5, "cm"),
    legend.margin = margin(0,0,0,0, "cm"),
    ## Strips for facet_wrap
    strip.text=element_text(size=8,family="sans",face="bold"),
    #strip.background=element_rect(fill="#f0f0f0")
    strip.background=element_blank()
    )



AAAS_palette <- c("blue1"="#3B4992FF","red1"="#EE0000FF","green1"="#008B45FF",
                  "purple1"="#631879FF","teal1"="#008280FF","red2"="#BB0021FF",
                  "purple2"="#5F559BFF","purple3"="#A20056FF",
                  "grey1"="#808180FF","black"="#1B1919FF")


theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){
  
  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space
  
  
  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))
  
  #Add axis
  theme.list <- 
    list(
      theme_void(), #Empty the current theme
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
      annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
      annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
      xlim(0 - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(0, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)
  
  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){
    
    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))
    
    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))
    
    #Add ticks to geom line for x axis
    #theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
    #                                         data = xtick, size = size, 
    #                                         color = color)
    
    #Add labels to the x-ticks
    #theme.list[[nlist + 4*k-2]] <- annotate("text", 
    #                                        x = ticks_x[k], 
    #                                        y = ygeo - 2.5*epsilon,
    #                                        size = textsize,
    #                                        label = paste(ticks_x[k]))
    
    
    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)
    
    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                          x = xgeo - 2.5*epsilon, 
                                          y = ticks_y[k],
                                          size = textsize,
                                          label = paste(ticks_y[k]))
  }
  
  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}

plot_nursing_home <- function(chain, obs_dat, plot_pars=c("viral_peak","R0","t0"), loc="NH 1"){
  
  
  ## Plot trajectories
  p_trajectories <- ggplot(posterior_dat) +
    geom_line(aes(x=t,y=40-prob_infection/0.004,group=sampno,col="Posterior draw"),size=0.1) +
    geom_line(data=best_prob_dat,aes(x=t,y=40-prob_infection/0.004,col="MAP"),size=0.5) +
    geom_violin(data=obs_dat,aes(x=t, y=ct,group=t),
                width=4,scale="width",fill="grey70",alpha=0.25,col="grey70",
                draw_quantiles=c(0.025,0.5,0.975)) +
    geom_jitter(data=obs_dat,aes(x=t, y=ct),size=0.25,height=0,width=1.25) +
    scale_y_continuous(trans="reverse",
                       sec.axis=sec_axis(~(40*0.004) - .*0.004,name="Per capita incidence",
                                         breaks=seq(0,0.12,by=0.03))) +
    xlab("Time") +
    ylab("Probability of infection") +
    theme_classic() + +
    scale_color_manual(values=c("Posterior draw"="gray50","MAP"="green","Ground truth"="blue")) +
    scale_linetype_manual(values=c("Sample date"="dashed")) +
    guides(col=guide_legend(title=NULL),linetype=guide_legend(title=NULL)) +
    theme(legend.position="bottom") +
    facet_wrap(~t, nrow=1)
  
  ggplot() +
    #geom_vline(data=betas_all_models_comb %>% dplyr::select(loc, t) %>% distinct() %>%
    #             rename(location=loc),
    #           aes(xintercept=t),col=AAAS_palette["grey1"],linetype="dashed",size=0.25) +
   
    
    geom_line(data=dat_use %>% group_by(location,mean_week_date) %>% summarize(y=median(ct, na.rm=TRUE)),
              aes(x=mean_week_date,y=y),col=AAAS_palette["grey"],linetype="longdash",size=0.25) +
    
    geom_ribbon(data=quants_inc,
                aes(ymin=40 - lower/0.004,ymax=40 - upper/0.004,x=time),
                alpha=0.1,fill=AAAS_palette["blue1"]) +
    geom_line(data=best_inc_traj, 
              aes(x=time,y=40 - inc/0.004),
              col=AAAS_palette["blue1"]) +
    
    scale_x_date(limits=range(gr_quants$time), breaks="7 days") +
    scale_y_continuous(trans="reverse",
                       sec.axis=sec_axis(~(40*0.004) - .*0.004,name="Per capita incidence",
                                         breaks=seq(0,0.12,by=0.03))) +
    ylab("Ct value") +
    xlab("Date") +
    export_theme
  
  ## Plot 
  
  ## Plot densities
  p_densities <- chain[,c("sampno",unique(parTab[which(parTab$fixed == 0),"names"]),"chain")] %>%
    mutate(chain = as.factor(chain)) %>%
    pivot_longer(-c(sampno,chain)) %>%
    ggplot() +
    geom_density(aes(x=value,fill=chain),alpha=0.25) +
    ggsci::scale_fill_aaas() +
    facet_wrap(t~name,scales="free",nrow=1) +
    export_theme
  
}


