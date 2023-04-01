library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(mgcv)
library(patchwork)
library(data.table)
source('seasonality_simulation/DGP_simulation_utils.R')
source('seasonality_dashboard/utils.R')

beta_int_seq <- seq(-0.4,0.4,by=0.05)
effect_simulation <- lapply(beta_int_seq, function(beta_int){
  print(beta_int)
  lapply(1:20,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(0.20,beta_int),phi=pi/4,a=a,b=b,mu_lag=1,sigma_lag=1)
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

effect_true_onset_plot <- pivot_longer(effect_simulation,cols=c('beta_binary','beta_qt')) %>%
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
                          geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
                          #geom_point(size=2) +
                          geom_abline(intercept=0,slope=1) +
                          #ggtitle('True vs fitted effect size for true disease onset data') +
                          #scale_x_continuous(breaks=round(beta_int_seq,2)) +
                          #scale_y_continuous(breaks=round(beta_int_seq,2)) +
                          scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
                          xlab(expression(beta[g])) +
                          ylab(expression(hat(beta)[pheno])) +
                          theme_bw()
ggsave('figures/simulation_presentation/effect_true_onset.png',effect_true_onset_plot,width=7,height=4)


pval_true_onset_plot <- pivot_longer(effect_simulation,cols=c('pval_binary','pval_qt')) %>%
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
                        geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
                        geom_abline(intercept=-log10(5e-8),slope=0,linetype='dashed') +
                        #ggtitle('P-values for true disease onset data') +
                        scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
                        xlab(expression(beta[g])) +
                        ylab(expression(-log[10](p))) +
                        theme_bw()
ggsave('figures/simulation_presentation/pval_true_onset.png',pval_true_onset_plot,width=7,height=4)

true_onset_plot <- (effect_true_onset_plot + theme(legend.position='none')) + pval_true_onset_plot
ggsave('figures/simulation_presentation/true_onset.png',true_onset_plot,width=8,height=4)


beta_baseline_seq <- seq(-0.5,0.5,by=0.1)
baseline_effect_simulation <- lapply(beta_baseline_seq, function(beta_base){
  print(beta_base)
  lapply(1:20,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(beta_base,-0.2),phi=pi/4,a=a,b=b,mu_lag=1,sigma_lag=1)
    tibble(sim_nr=i,
           beta_base=beta_base,
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

effect_true_onset_baseline_plot <- pivot_longer(baseline_effect_simulation,cols=c('beta_binary','beta_qt')) %>%
                                  mutate(name=gsub('beta_','',name)) %>%
                                  group_by(name,beta_base) %>%
                                  summarise(avg_effect=mean(value),
                                            lower_effect=quantile(value,0.025),
                                            upper_effect=quantile(value,0.975)) %>%
                                  ggplot(aes(x=beta_base,
                                             y=avg_effect,
                                             ymin=lower_effect,
                                             ymax=upper_effect,
                                             col=factor(name))) +
                                  geom_pointrange(position = position_jitterdodge(dodge.width = 0.005, jitter.width = 0)) +
                                  #ggtitle('True vs fitted effect size for true disease onset data') +
                                  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
                                  xlab(expression(beta[baseline])) +
                                  ylab(expression(hat(beta)[pheno])) +
                                  theme_bw()
ggsave('figures/simulation_presentation/effect_true_onset_baseline.png',effect_true_onset_baseline_plot,width=7,height=4)


pval_true_onset_baseline_plot <- pivot_longer(baseline_effect_simulation,cols=c('pval_binary','pval_qt')) %>%
                                  mutate(name=gsub('pval_','',name),
                                         value=-log10(value)) %>%
                                  group_by(name,beta_base) %>%
                                  summarise(avg_effect=mean(value),
                                            lower_effect=quantile(value,0.025),
                                            upper_effect=quantile(value,0.975)) %>%
                                  ggplot(aes(x=beta_base,
                                             y=avg_effect,
                                             ymin=lower_effect,
                                             ymax=upper_effect,
                                             col=factor(name))) +
                                  geom_pointrange(position = position_jitterdodge(dodge.width = 0.005, jitter.width = 0)) +
                                  geom_abline(intercept=-log10(5e-8),slope=0) +
                                  #ggtitle('P-values for true disease onset data') +
                                  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
                                  xlab(expression(beta[baseline])) +
                                  ylab(expression(-log[10](p))) +
                                  theme_bw()
ggsave('figures/simulation_presentation/pval_true_onset_baseline.png',pval_true_onset_baseline_plot,width=7,height=4)

true_onset_baseline_plot <- (effect_true_onset_baseline_plot + theme(legend.position='none')) + pval_true_onset_baseline_plot
ggsave('figures/simulation_presentation/true_onset_baseline.png',true_onset_baseline_plot,width=8,height=4)

mu_lag_seq <- seq(0,6,by=0.5)
a <- 1;b <- 13
lag_simulation <- lapply(mu_lag_seq, function(mu_lag){
  print(mu_lag)
  lapply(1:20,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(0.2,-0.3),phi=pi/4,mod_type='additive',a=a,b=b,mu_lag=mu_lag,sigma_lag=sqrt(mu_lag)+1)
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


effect_lag_plot <- pivot_longer(lag_simulation,cols=c('beta_binary_pb','beta_qt_pb')) %>%
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
                  #ggtitle('Effect size for perturbed disease onset data') +
                  #scale_x_continuous(breaks=round(beta_int_seq,2)) +
                  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
                  xlab('Mean diagnosis lag (in months)') +
                  ylab(expression(hat(beta)[pheno])) +
                  theme_bw()
ggsave('figures/simulation_presentation/effect_pb_lag.png',effect_lag_plot,width=7,height=4)

pval_lag_plot <- pivot_longer(lag_simulation,cols=c('pval_binary_pb','pval_qt_pb')) %>%
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
                  #ggtitle('P-value of associations on perturbed disease onset data') +
                  scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
                  xlab('Mean diagnosis lag (in months)') +
                  ylab(expression(-log[10](p))) +
                  theme_bw()
ggsave('figures/simulation_presentation/pval_pb_lag.png',pval_lag_plot,width=7,height=4)

lag_plot <- (effect_lag_plot + theme(legend.position='none')) + pval_lag_plot
ggsave('figures/simulation_presentation/lag_plot.png',lag_plot,width=8,height=4)



sim_dat <- simulate_seasonal_dist(n=40000,MAF=0.02,beta=c(0.2,-0.15),phi=pi/4,mod_type='additive',a=a,b=b)
mod_obj <- seasonal_mod(t=sim_dat$EVENT_MONTH_DEC,x_g=sim_dat$gt_mod,a=1,b=13)


