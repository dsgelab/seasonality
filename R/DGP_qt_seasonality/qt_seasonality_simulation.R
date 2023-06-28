library(ggplot2)
library(dplyr)
library(mgcv)


dgp_mu <- function(beta,phi,x_g,t,P){
  beta['g']*x_g + (beta['a'] + beta['ag']*x_g)*sin(t*2*pi/P + phi + beta['pg']*x_g)
}

get_mu_grid <- function(beta,phi,mod_type,P){
  t_seq <- seq(0,P,by=0.01)
  gt_mod <- 0:2
  gt_lab <- 0:2
  if(mod_type=='dominant'){
    gt_lab <- gt_mod[1:2]
    gt_mod <- gt_mod[1:2]
  }else if(mod_type=='recessive'){
    gt_lab <- gt_mod[c(1,3)]
    gt_mod <- gt_mod[1:2]
  }
  mu_grid <- lapply(1:length(gt_mod),function(i){
                tibble(t=t_seq,
                       gt_mod=gt_mod[i],
                       gt_lab=gt_lab[i],
                       mu=dgp_mu(beta=beta,phi=phi,x_g=gt_mod[i],t=t_seq,P=P))
              }) %>% bind_rows()
  return(mu_grid)
}

simulate_seasonal_qt <- function(n,beta,phi,sigma,MAF,mod_type,P){
  t <- P*runif(n) # assume time of measurements are uniformly spread over the interval [0,P)
  x_g <- rbinom(n,size=2,prob=MAF) # simulate genotype in Harvey-Weinberg equilibrium
  if(mod_type=='dominant'){
    gt_lab <- as.integer(x_g>0)
    gt_mod <- gt_lab
  }else if(mod_type=='recessive'){
    gt_lab <- 2*as.integer(x_g>1)
    gt_mod <- as.integer(x_g>1)
  }else{
    gt_mod <- x_g
    gt_lab <- x_g
  }
  mu <- dgp_mu(beta,phi,gt_mod,t,P)
  y <- rnorm(n=n,mean=mu,sd=sigma)
  sim_dat <- tibble(t=t,gt_mod=gt_mod,gt_lab=gt_lab,y=y)
  #adjust genotype labels based on model type
  return(sim_dat)
}
run_seasonal_qt_gam <- function(dat,P){
  f_seasonal <- paste('y ~ s(t, k=6, bs="cp")',sep='+')
  mod_seasonal <- gam(as.formula(f_seasonal), family=gaussian(), data=dat, knots=list(t = c(0, P)), scale=-1)
  return(mod_seasonal)
}

add_seasonal_values <- function(dat,mod_gam,P){
  t_seq <- seq(0,P,by=0.01)
  seasonal_pred <- predict(mod_gam,newdata=data.frame(t=t_seq),type='terms',se.fit=T)
  seasonal_pred_dat <- tibble(t=t_seq,
                              est=seasonal_pred$fit[,'s(t)'],
                              lower=est - 1.96*seasonal_pred$se.fit[,'s(t)'],
                              upper=est + 1.96*seasonal_pred$se.fit[,'s(t)'])
    out_dat <- mutate(dat,
                      seasonal_val=approx(x=seasonal_pred_dat$t,
                                          y=seasonal_pred_dat$est,
                                          xout=dat$t)$y,
                      seasonal_binary=seasonal_val>0)
    return(out_dat)
}

run_seasonal_amplitude_mod_qt <- function(dat){
  high_season_mod <- lm(y~gt_mod,data=filter(dat,seasonal_binary>0))
  low_season_mod <- lm(y~gt_mod,data=filter(dat,seasonal_binary<=0))

  seasonal_binary_summary <- tibble(effect=coef(high_season_mod)['gt_mod']-coef(low_season_mod)['gt_mod'],
                             se=sqrt(summary(high_season_mod)$coefficients['gt_mod',2]^2+summary(low_season_mod)$coefficients['gt_mod',2]^2), #add squares of se's to get effect se
                             lower=effect-1.96*se,
                             upper=effect+1.96*se,
                             wald=(effect/se)^2,
                             pval=pchisq(wald,df=1,ncp=0,lower.tail=F))

  return(seasonal_binary_summary)
}

# Plotting
plot_mu_grid <- function(mu_grid,mod_type,P){
  cols <- c("#BC3C29FF","#0072B5FF","#E18727FF")
  if(mod_type=='dominant'){
    cols <- cols[1:2]
  }else if(mod_type=='recessive'){
    cols <- cols[c(1,3)]
  }
  ggplot(data=mu_grid) +
    geom_line(aes(t,mu,col=factor(gt_lab))) +
    scale_color_manual(values=cols,name='Genotype') +
    scale_x_continuous(breaks=seq(0,P)) +
    ylab('Expected value') +
    xlab('Month number') +
    theme_classic()
}

plot_sim_dat <- function(sim_dat,mu_grid,mod_type,P){
  plot_mu_grid(mu_grid=mu_grid,mod_type=mod_type,P=P) +
  geom_point(data=sim_dat,aes(t,y,col=factor(gt_lab)),alpha=0.1) +
  ylab('y')
}

plot_seasonality_qt <- function(mod,P){
  t_seq <- seq(0,P,by=0.01)
  seasonal_pred <- predict(mod,newdata=data.frame(t=t_seq),type='terms',se.fit=T)
  seasonal_pred_dat <- tibble(t=t_seq,
                              est=seasonal_pred$fit[,'s(t)'],
                              lower=est - 1.96*seasonal_pred$se.fit[,'s(t)'],
                              upper=est + 1.96*seasonal_pred$se.fit[,'s(t)'])
  ggplot(seasonal_pred_dat) +
    geom_line(aes(x=t,
                  y=est)) +
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper),alpha=0.3,fill='blue') +
    geom_hline(yintercept=0,linetype='dashed',color='red') +
    scale_x_continuous(breaks=seq_len(12)) +
    xlab('Month number') +
    ylab('Smooth term') +
    ggtitle('Seasonal component') +
    theme_bw() +
    theme(axis.title.x=element_text(size=14),
          axis.text.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          axis.text.y=element_text(size=14))
}

# Assumptions about the data-generating mechanism
input <- list(n=10000,
              MAF=0.1,
              P=12, # period
              beta=c('g'=0,  # additive genetic effect
                     'a'=0.2,  # baseline seasonal/circadian amplitude
                     'ag'=0.1,   # genetic-season amplitude interaction
                     'pg'=0), # genetic-season phase interaction
              phi=1, # baseline phase
              sigma=0.05, #residual standard deviation
              mod_type='additive')

mu_grid <- get_mu_grid(beta=input$beta,phi=input$phi,mod_type=input$mod_type,P=12)
#Visualize theoretical mean seasonality curves
plot_mu_grid(mu_grid,mod_type=input$mod_type,P=12)

sim_dat <- simulate_seasonal_qt(n=input$n,beta=input$beta,phi=input$phi,sigma=input$sigma,MAF=input$MAF,mod_type=input$mod_type,P=12)
#Scatter plot of the simulated data with the theoretical mean curve on top
plot_sim_dat(sim_dat=sim_dat,mu_grid=mu_grid,mod_type=input$mod_type,P=12)

#Run GAM to extract baseline seasonality and plot
mod_gam <- run_seasonal_qt_gam(sim_dat,P=12)
# Visualize the baseline seasonality in the simulated data
plot_seasonality_qt(mod_gam,P=12)

#add baseline seasonal information for seasonality GWAS (comparing high season vs. low season)
sim_dat <- add_seasonal_values(sim_dat,mod_gam,P=12)
run_seasonal_amplitude_mod_qt(dat=sim_dat)

#demonstrate seasonal splitting
ggplot(sim_dat) +
geom_point(aes(t,y,col=factor(seasonal_binary)),alpha=0.5) +
scale_color_manual(values=c("#BC3C29FF","#0072B5FF"),name='High season') +
theme_classic()

