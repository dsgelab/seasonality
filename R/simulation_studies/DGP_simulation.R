library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(mgcv)
library(patchwork)
library(data.table)
source('R/simulation_studies/DGP_simulation_utils.R')
source('seasonality_dashboard/utils.R')


tibble(t=c(0,0.999,1,12.999,13,20),dens=c(0,0,1/12,1/12,0,0)) %>%
ggplot(aes(t,dens)) +
geom_line(size=1,color='#E18727FF') +
scale_y_continuous(limits=c(0,0.15)) +
xlab('t') +
ylab('Density') +
theme_classic()

ggsave('figures/simulation_presentation/f_no_seasonality.png',width=5,height=3)


tibble(t=c(0,0.999,1,12.999,13,20),surv=c(0,0,1,0,0,0)) %>%
ggplot(aes(t,surv)) +
geom_line(size=1,color='#E18727FF') +
xlab('t') +
ylab('Survival') +
theme_classic()

ggsave('figures/simulation_presentation/S_no_seasonality.png',width=5,height=3)
a <- 1;b <- 13
dist_grid_seasonality <- get_dist_grid(beta=c(0.5,0),phi=pi/2,a=a,b=b) %>%
                         bind_rows() %>%
                         filter(gt==0)


seasonality_d <- ggplot(dist_grid_seasonality,aes(t,dens)) +
                 geom_line(size=1,color='#E18727FF') +
                 xlab('t') +
                 ylab('Density') +
                 theme_classic()

seasonality_S <- ggplot(dist_grid_seasonality,aes(t,surv)) +
                 geom_line(size=1,color='#E18727FF') +
                 xlab('t') +
                 ylab('Survival') +
                 theme_classic()

comb_seasonality <- seasonality_d + seasonality_S

ggsave('figures/simulation_presentation/comb_seasonality.png',comb_seasonality,width=7,height=3)


dist_grid_seasonality_int <- get_dist_grid(beta=c(0.5,-0.3),phi=pi/2,a=a,b=b) %>%
                             bind_rows()

seasonality_d_int <- ggplot(dist_grid_seasonality_int,aes(t,dens,col=as.factor(gt))) +
                      geom_line(size=1) +
                      scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='Genotype') +
                      xlab('t') +
                      ylab('Density') +
                      theme_classic()

seasonality_S_int <- ggplot(dist_grid_seasonality_int,aes(t,surv,col=as.factor(gt))) +
                    geom_line(size=1) +
                    scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='Genotype') +
                    xlab('t') +
                    ylab('Survival') +
                    theme_classic() +
                    theme(legend.position = 'none')



comb_seasonality_int <- seasonality_d_int + seasonality_S_int

ggsave('figures/simulation_presentation/comb_seasonality_int.png',comb_seasonality_int,width=10,height=3)


sim_dat <- simulate_seasonal_dist(n=40000,MAF=0.1,beta=c(0.25,-0.25),phi=pi/2,a=a,b=b)

surv_obj_sim <- Surv(time=sim_dat$EVENT_MONTH_DEC,event=rep(1,nrow(sim_dat)))
survfit_sim <- survfit(surv_obj_sim ~ sim_dat$gt)

KM_gg(survfit_sim,'True disease onset')
ggsave('figures/simulation_presentation/KM_sim.png',width=5,height=3)

#+ (KM_gg(survfit_sim_pb,title='True disease onset with perturbation') + theme(legend.position='none'))

monthly_counts <- mutate(sim_dat,EVENT_MONTH=floor(EVENT_MONTH_DEC)) %>%
                  count(EVENT_MONTH,name = 'count') %>%
                  mutate(ENDPOINT_LAB='True disease onset')
ggplot(monthly_counts,aes(EVENT_MONTH,count)) +
geom_col(position='dodge',fill="#0072B5FF",width=0.7) +
scale_x_continuous(breaks=seq(a,b),expand=c(0,0.2)) +
scale_y_continuous(expand=c(0,0)) +
xlab('Month') +
ylab('Count') +
theme_bw()
ggsave('figures/simulation_presentation/monthly_counts.png',width=5,height=3)

mod_counts <- run_seasonality_gam_sim(monthly_counts)
seasonality_plot(mod_dat,list(seasonal=mod_counts))
ggsave('figures/simulation_presentation/seasonality_spline.png',width=5,height=3)


beta_int_seq <- seq(-0.3,0.3,by=0.05)
effect_simulation <- lapply(beta_int_seq, function(beta_int){
  print(beta_int)
  lapply(1:10,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(0.25,beta_int),phi=pi/2,a=a,b=b,mu_lag=1,sigma_lag=2)
    tibble(sim_nr=i,
           beta_int=beta_int,
           beta_binary=summary(sim$binary_mod)$coefficients[2,1],
           se_binary=summary(sim$binary_mod)$coefficients[2,2],
           pval_binary=summary(sim$binary_mod)$coefficients[2,4],
           beta_qt=summary(sim$qt_mod)$coefficients[2,1],
           se_qt=summary(sim$qt_mod)$coefficients[2,2],
           pval_qt=summary(sim$qt_mod)$coefficients[2,4],
           beta_binary_pb=summary(sim$binary_mod_pb)$coefficients[2,1],
           se_binary_pb=summary(sim$binary_mod_pb)$coefficients[2,2],
           pval_binary_pb=summary(sim$binary_mod_pb)$coefficients[2,4],
           beta_qt_pb=summary(sim$qt_mod_pb)$coefficients[2,1],
           se_qt_pb=summary(sim$qt_mod_pb)$coefficients[2,2],
           pval_qt_pb=summary(sim$qt_mod_pb)$coefficients[2,4])
  }) %>% bind_rows()
}) %>% bind_rows()

pivot_longer(effect_simulation,cols=c('beta_binary','beta_qt')) %>%
  mutate(name=gsub('beta_','',name)) %>%
  group_by(name,beta_int) %>%
  summarise(avg_effect=mean(value),
            lower_effect=quantile(value,0.025),
            upper_effect=quantile(value,0.975)) %>%
  ggplot(aes(x=beta_int,
             y=avg_effect,
             ymin=lower_effect,
             ymax=upper_effect,
             col=factor(name))) +
  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.005, jitter.width = 0)) +
  #geom_point(size=2) +
  geom_abline(intercept=0,slope=1) +
  ggtitle('True vs fitted effect size for true disease onset data') +
  #scale_x_continuous(breaks=round(beta_int_seq,2)) +
  #scale_y_continuous(breaks=round(beta_int_seq,2)) +
  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
  xlab(expression(beta[g])) +
  ylab(expression(hat(beta)[pheno])) +
  theme_bw()
ggsave('figures/simulation_presentation/effect_true_onset.png',width=7,height=4)


pivot_longer(effect_simulation,cols=c('pval_binary','pval_qt')) %>%
  mutate(name=gsub('pval_','',name),
         value=-log10(value)) %>%
  group_by(name,beta_int) %>%
  summarise(avg_effect=mean(value),
            lower_effect=quantile(value,0.025),
            upper_effect=quantile(value,0.975)) %>%
  ggplot(aes(x=beta_int,
             y=avg_effect,
             ymin=lower_effect,
             ymax=upper_effect,
             col=factor(name))) +
  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.005, jitter.width = 0)) +
  geom_abline(intercept=-log10(5e-8),slope=0) +
  ggtitle('P-values as a function of true effect size for true disease onset data') +
  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
  xlab(expression(beta[g])) +
  ylab(expression(-log[10](p))) +
  theme_bw()
ggsave('figures/simulation_presentation/pval_true_onset.png',width=7,height=4)


pivot_longer(effect_simulation,cols=c('beta_binary','beta_binary_pb')) %>%
  mutate(name=gsub('beta_','',name)) %>%
  group_by(name,beta_int) %>%
  summarise(avg_effect=mean(value),
            lower_effect=quantile(value,0.025),
            upper_effect=quantile(value,0.975)) %>%
  ggplot(aes(x=beta_int,
             y=avg_effect,
             ymin=lower_effect,
             ymax=upper_effect,
             col=factor(name))) +
  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.005, jitter.width = 0)) +
  #geom_point(size=2) +
  geom_abline(intercept=0,slope=1) +
  ggtitle('True vs fitted effect size for true and perturbed disease onset data') +
  #scale_x_continuous(breaks=round(beta_int_seq,2)) +
  #scale_y_continuous(breaks=round(beta_int_seq,2)) +
  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
  xlab(expression(beta[g])) +
  ylab(expression(hat(beta)[pheno])) +
  theme_bw()
ggsave('figures/simulation_presentation/effect_true_pb_onset.png',width=7,height=4)


pivot_longer(effect_simulation,cols=c('pval_binary','pval_binary_pb')) %>%
  mutate(name=gsub('pval_','',name),
         value=-log10(value)) %>%
  group_by(name,beta_int) %>%
  summarise(avg_effect=mean(value),
            lower_effect=quantile(value,0.025),
            upper_effect=quantile(value,0.975)) %>%
  ggplot(aes(x=beta_int,
             y=avg_effect,
             ymin=lower_effect,
             ymax=upper_effect,
             col=factor(name))) +
  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.005, jitter.width = 0)) +
  geom_abline(intercept=-log10(5e-8),slope=0) +
  ggtitle('P-values for true and perturbed disease onset data') +
  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
  xlab(expression(beta[g])) +
  ylab(expression(-log[10](p))) +
  theme_bw()
ggsave('figures/simulation_presentation/pval_true_pb_onset.png',width=7,height=4)

mu_lag_seq <- seq(0,6,by=0.5)
lag_simulation <- lapply(mu_lag_seq, function(mu_lag){
  print(mu_lag)
  lapply(1:10,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(0.25,-0.5),phi=pi/2,a=a,b=b,mu_lag=mu_lag,sigma_lag=2+mu_lag/3)
    tibble(sim_nr=i,
           mu_lag=mu_lag,
           beta_binary=summary(sim$binary_mod)$coefficients[2,1],
           se_binary=summary(sim$binary_mod)$coefficients[2,2],
           pval_binary=summary(sim$binary_mod)$coefficients[2,4],
           beta_qt=summary(sim$qt_mod)$coefficients[2,1],
           se_qt=summary(sim$qt_mod)$coefficients[2,2],
           pval_qt=summary(sim$qt_mod)$coefficients[2,4],
           beta_binary_pb=summary(sim$binary_mod_pb)$coefficients[2,1],
           se_binary_pb=summary(sim$binary_mod_pb)$coefficients[2,2],
           pval_binary_pb=summary(sim$binary_mod_pb)$coefficients[2,4],
           beta_qt_pb=summary(sim$qt_mod_pb)$coefficients[2,1],
           se_qt_pb=summary(sim$qt_mod_pb)$coefficients[2,2],
           pval_qt_pb=summary(sim$qt_mod_pb)$coefficients[2,4])
  }) %>% bind_rows()
}) %>% bind_rows()


pivot_longer(lag_simulation,cols=c('beta_binary_pb','beta_qt_pb')) %>%
  mutate(name=gsub('beta_','',name)) %>%
  group_by(name,mu_lag) %>%
  summarise(avg_effect=mean(value),
            lower_effect=quantile(value,0.025),
            upper_effect=quantile(value,0.975)) %>%
  ggplot(aes(x=mu_lag,
             y=avg_effect,
             ymin=lower_effect,
             ymax=upper_effect,
             col=factor(name))) +
  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0)) +
  geom_abline(intercept=0,slope=0,linetype='longdash') +
  ggtitle('Effect size for perturbed disease onset data') +
  #scale_x_continuous(breaks=round(beta_int_seq,2)) +
  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
  xlab('Mean diagnosis lag (in months)') +
  ylab(expression(hat(beta)[pheno])) +
  theme_bw()
ggsave('figures/simulation_presentation/effect_pb_lag.png',width=7,height=4)

pivot_longer(lag_simulation,cols=c('pval_binary_pb','pval_qt_pb')) %>%
  mutate(name=gsub('beta_','',name),
         value=-log10(value)) %>%
  group_by(name,mu_lag) %>%
  summarise(avg_effect=mean(value),
            lower_effect=quantile(value,0.025),
            upper_effect=quantile(value,0.975)) %>%
  ggplot(aes(x=mu_lag,
             y=avg_effect,
             ymin=lower_effect,
             ymax=upper_effect,
             col=factor(name))) +
  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0)) +
  geom_abline(intercept=-log10(5e-8),slope=0,linetype='longdash') +
  ggtitle('P-value of associations on perturbed disease onset data') +
  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
  xlab('Mean diagnosis lag (in months)') +
  ylab(expression(-log[10](p))) +
  theme_bw()
ggsave('figures/simulation_presentation/pval_pb_lag.png',width=7,height=4)

eval_g <- function(theta){
  c(-1+1e-6 + theta[1] + 2*theta[2],-1+1e-6 - theta[1] - 2*theta[2])
}

eval_jac_g <- function(theta) {
  return(rbind(c(1,2,0),c(-1,-2,0)))
}

local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 1000,
              "local_opts" = local_opts,
              "print_level" = 0 )
res <- nloptr(x0 = c(0,0,0),eval_f = loglik_seasonal,lb = c(-1+1e-6,-0.5+1e-6,0), ub = c(1,0.5,pi)-1e-6,eval_g_ineq = eval_g,opts = opts)

