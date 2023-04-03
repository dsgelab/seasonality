library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(mgcv)
library(patchwork)
library(data.table)
source('seasonality_dashboard/DGP_simulation_utils.R')
source('seasonality_dashboard/utils.R')
a <- 1;b <- 13
beta_int_seq <- seq(-0.4,0.4,by=0.05)
genetic_effect_simulation <- lapply(beta_int_seq, function(beta_int){
  print(beta_int)
  lapply(1:20,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(0.20,beta_int),phi=pi/4,a=a,b=b,mu_lag=1,sigma_lag=1)
    mutate(sim$mod_summary,sim_nr=i,beta_int=beta_int) %>%
    select(sim_nr,beta_int,everything())
  }) %>% bind_rows()
}) %>% bind_rows()

write.table(genetic_effect_simulation,'results/simulations/genetic_effect_simulation.tsv',sep='\t',row.names=F,quote=F)

effect_true_onset_plot <- filter(genetic_effect_simulation,!perturbed) %>%
                          group_by(mod,beta_int) %>%
                          summarise(avg_effect=mean(estimate),
                                    lower_effect=quantile(estimate,0.025),
                                    upper_effect=quantile(estimate,0.975)) %>%
                          ggplot(aes(x=beta_int,
                                     y=avg_effect,
                                     ymin=lower_effect,
                                     ymax=upper_effect,
                                     col=factor(mod))) +
                          geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
                          geom_abline(intercept=0,slope=1) +
                          scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                          xlab(expression(beta[g])) +
                          ylab(expression(hat(beta)[pheno])) +
                          theme_bw()


pval_true_onset_plot <- filter(genetic_effect_simulation,!perturbed) %>%
                        group_by(mod,beta_int) %>%
                        mutate(minus_log10_p=-log10(p.value)) %>%
                        summarise(avg_effect=mean(minus_log10_p),
                                  lower_effect=quantile(minus_log10_p,0.025),
                                  upper_effect=quantile(minus_log10_p,0.975)) %>%
                                              ggplot(aes(x=beta_int,
                                                         y=avg_effect,
                                                         ymin=lower_effect,
                                                         ymax=upper_effect,
                                                         col=factor(mod))) +
                        geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
                        geom_abline(intercept=-log10(5e-8),slope=0,linetype='dashed') +
                        scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                        xlab(expression(beta[g])) +
                        ylab(expression(-log[10](p))) +
                        theme_bw()

true_onset_plot <- (effect_true_onset_plot + theme(legend.position='none')) + pval_true_onset_plot
ggsave('figures/power_simulation/effect_pval_true_onset.png',true_onset_plot,width=8,height=4)


beta_baseline_seq <- seq(0.05,0.5,by=0.05)
baseline_effect_simulation <- lapply(beta_baseline_seq, function(beta_base){
  print(beta_base)
  lapply(1:20,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(beta_base,-0.2),phi=pi/4,a=a,b=b,mu_lag=1,sigma_lag=1)
    mutate(sim$mod_summary,sim_nr=i,beta_base=beta_base) %>%
      select(sim_nr,beta_base,everything())
  }) %>% bind_rows()
}) %>% bind_rows()
write.table(baseline_effect_simulation,'results/simulations/baseline_effect_simulation.tsv',sep='\t',row.names=F,quote=F)

effect_true_onset_baseline_plot <- filter(baseline_effect_simulation,!perturbed) %>%
                                    group_by(mod,beta_base) %>%
                                    summarise(avg_effect=mean(estimate),
                                              lower_effect=quantile(estimate,0.025),
                                              upper_effect=quantile(estimate,0.975)) %>%
                                    ggplot(aes(x=beta_base,
                                               y=avg_effect,
                                               ymin=lower_effect,
                                               ymax=upper_effect,
                                               col=factor(mod))) +
                                    geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
                                    scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                                    xlab(expression(beta[b])) +
                                    ylab(expression(hat(beta)[pheno])) +
                                    theme_bw()


pval_true_onset_baseline_plot <- filter(baseline_effect_simulation,!perturbed) %>%
                                  group_by(mod,beta_base) %>%
                                  mutate(minus_log10_p=-log10(p.value)) %>%
                                  summarise(avg_effect=mean(minus_log10_p),
                                            lower_effect=quantile(minus_log10_p,0.025),
                                            upper_effect=quantile(minus_log10_p,0.975)) %>%
                                  ggplot(aes(x=beta_base,
                                             y=avg_effect,
                                             ymin=lower_effect,
                                             ymax=upper_effect,
                                             col=factor(mod))) +
                                  geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
                                  geom_abline(intercept=-log10(5e-8),slope=0,linetype='dashed') +
                                  scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                                  xlab(expression(beta[b])) +
                                  ylab(expression(-log[10](p))) +
                                  theme_bw()

true_onset_baseline_plot <- (effect_true_onset_baseline_plot + theme(legend.position='none')) + pval_true_onset_baseline_plot
ggsave('figures/power_simulation/effect_pval_true_onset_baseline.png',true_onset_baseline_plot,width=8,height=4)

mu_lag_seq <- seq(0,6,by=0.5)
a <- 1;b <- 13
lag_simulation <- lapply(mu_lag_seq, function(mu_lag){
  print(mu_lag)
  lapply(1:20,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(0.2,-0.3),phi=pi/4,mod_type='additive',a=a,b=b,mu_lag=mu_lag,sigma_lag=sqrt(mu_lag)+1)
    mutate(sim$mod_summary,sim_nr=i,mu_lag=mu_lag) %>%
      select(sim_nr,mu_lag,everything())
  }) %>% bind_rows()
}) %>% bind_rows()
write.table(lag_simulation,'results/simulations/lag_simulation.tsv',sep='\t',row.names=F,quote=F)


effect_lag_plot <- filter(lag_simulation,perturbed) %>%
                    group_by(mod,mu_lag) %>%
                    summarise(avg_effect=mean(estimate),
                              lower_effect=quantile(estimate,0.025),
                              upper_effect=quantile(estimate,0.975)) %>%
                                    ggplot(aes(x=mu_lag,
                                               y=avg_effect,
                                               ymin=lower_effect,
                                               ymax=upper_effect,
                                               col=factor(mod))) +
                  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0)) +
                  geom_abline(intercept=0,slope=0,linetype='longdash') +
                  #ggtitle('Effect size for perturbed disease onset data') +
                  #scale_x_continuous(breaks=round(beta_int_seq,2)) +
                  scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                  xlab('Mean diagnosis lag (in months)') +
                  ylab(expression(hat(beta)[pheno])) +
                  theme_bw()

pval_lag_plot <- filter(lag_simulation,perturbed) %>%
                group_by(mod,mu_lag) %>%
                mutate(minus_log10_p=-log10(p.value)) %>%
                summarise(avg_effect=mean(minus_log10_p),
                          lower_effect=quantile(minus_log10_p,0.025),
                          upper_effect=quantile(minus_log10_p,0.975)) %>%
                ggplot(aes(x=mu_lag,
                           y=avg_effect,
                           ymin=lower_effect,
                           ymax=upper_effect,
                           col=factor(mod))) +
                  geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0)) +
                  geom_abline(intercept=-log10(5e-8),slope=0,linetype='longdash') +
                  #ggtitle('P-value of associations on perturbed disease onset data') +
                  scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                  xlab('Mean diagnosis lag (in months)') +
                  ylab(expression(-log[10](p))) +
                  theme_bw()

lag_plot <- (effect_lag_plot + theme(legend.position='none')) + pval_lag_plot
ggsave('figures/power_simulation/effect_pval_lag_plot.png',lag_plot,width=8,height=4)



# sim_dat <- simulate_seasonal_dist(n=40000,MAF=0.02,beta=c(0.2,-0.15),phi=pi/4,mod_type='additive',a=a,b=b)
# mod_obj <- seasonal_mod(t=sim_dat$EVENT_MONTH_DEC,x_g=sim_dat$gt_mod,a=1,b=13)

## Power simulation
# vary effect size
effect_vec <- seq(0,0.5,by=0.05)
n_sim <- 100
effect_power_simulation <- lapply(effect_vec, function(beta_int){
  print(beta_int)
  lapply(1:100,function(i){
    sim <- sim_pipeline(n=40000,MAF=0.02,beta=c(-0.2,beta_int),phi=pi/4,mod_type='additive',a=a,b=b,mu_lag=0,sigma_lag=1)
    mutate(sim$mod_summary,sim_nr=i,beta_int=beta_int) %>%
    select(sim_nr,beta_int,everything())
  }) %>% bind_rows()
}) %>% bind_rows()
write.table(effect_power_simulation,'results/simulations/effect_power_simulation',sep='\t',row.names=F,quote=F)


effect_power <- filter(effect_power_simulation,!perturbed) %>%
                mutate(mod=factor(mod,levels=c('binary_mod','qt_mod','qt_raw_mod'),labels=c('Binary model','QT - normalized','QT - unnormalized'))) %>%
                group_by(mod,beta_int) %>%
                summarise(num_GWAS_signif=sum(p.value<5e-8),
                          power=num_GWAS_signif/n_sim,
                          power_se=sqrt(power*(1-power)/n_sim)) %>%
                ggplot() +
                geom_line(aes(beta_int,power,col=mod)) +
                geom_ribbon(aes(beta_int,ymin=power-1.96*power_se,ymax=power+1.96*power_se,fill=mod),alpha=0.1) +
                scale_x_continuous(expand=c(0,0.01)) +
                scale_y_continuous(expand=c(0,0.01)) +
                scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                scale_fill_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
                xlab(expression(beta[g])) +
                ylab('Power') +
                ggtitle('MAF=0.02, beta_b=-0.2, n=40000') +
                theme_classic() +
                theme(legend.position=c(0.8,0.3))
ggsave('figures/power_simulation/effect_power.png',effect_power,width=4,height=4)



#vary MAF
MAF_vec <- c(0.001,seq(0.01,0.05,by=0.01),seq(0.1,0.3,by=0.05))
n_sim <- 100
MAF_power_simulation <- lapply(MAF_vec, function(maf){
  print(maf)
  lapply(1:100,function(i){
    sim <- sim_pipeline(n=40000,MAF=maf,beta=c(-0.4,0.1),phi=pi/4,mod_type='additive',a=a,b=b,mu_lag=0,sigma_lag=1)
    mutate(sim$mod_summary,sim_nr=i,maf=maf) %>%
    select(sim_nr,maf,everything())
  }) %>% bind_rows()
}) %>% bind_rows()
write.table(MAF_power_simulation,'results/simulations/MAF_power_simulation.tsv',sep='\t',row.names=F,quote=F)

MAF_power <- filter(MAF_power_simulation,!perturbed) %>%
            mutate(mod=factor(mod,levels=c('binary_mod','qt_mod','qt_raw_mod'),labels=c('Binary model','QT - normalized','QT - unnormalized'))) %>%
            group_by(mod,maf) %>%
            summarise(num_GWAS_signif=sum(p.value<5e-8),
                      power=num_GWAS_signif/n_sim,
                      power_se=sqrt(power*(1-power)/n_sim)) %>%
            ggplot() +
            geom_line(aes(maf,power,col=mod)) +
            geom_ribbon(aes(maf,ymin=power-1.96*power_se,ymax=power+1.96*power_se,fill=mod),alpha=0.1) +
            scale_x_continuous(expand=c(0,0.01)) +
            scale_y_continuous(expand=c(0,0.01)) +
            scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
            scale_fill_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
            xlab('MAF') +
            ylab('Power') +
            ggtitle('beta_b=-0.4, beta_g=0.1, n=40000') +
            theme_classic() +
            theme(legend.position=c(0.8,0.3))
ggsave('figures/power_simulation/MAF_power.png',MAF_power,width=4,height=4)

#vary sample size
n_vec <- c(1000,5000,seq(10000,50000,by=10000))
n_sim <- 100
n_power_simulation <- lapply(n_vec, function(n){
  print(n)
  lapply(1:100,function(i){
    sim <- sim_pipeline(n=n,MAF=0.05,beta=c(-0.2,0.2),phi=pi/4,mod_type='additive',a=a,b=b,mu_lag=0,sigma_lag=1)
    mutate(sim$mod_summary,sim_nr=i,n=n) %>%
      select(sim_nr,n,everything())
  }) %>% bind_rows()
}) %>% bind_rows()
write.table(n_power_simulation,'results/simulations/n_power_simulation.tsv',sep='\t',row.names=F,quote=F)

n_power <- filter(n_power_simulation,!perturbed) %>%
            mutate(mod=factor(mod,levels=c('binary_mod','qt_mod','qt_raw_mod'),labels=c('Binary model','QT - normalized','QT - unnormalized'))) %>%
            group_by(mod,n) %>%
            summarise(num_GWAS_signif=sum(p.value<5e-8),
                      power=num_GWAS_signif/n_sim,
                      power_se=sqrt(power*(1-power)/n_sim)) %>%
            ggplot() +
            geom_line(aes(n,power,col=mod)) +
            geom_ribbon(aes(n,ymin=power-1.96*power_se,ymax=power+1.96*power_se,fill=mod),alpha=0.1) +
            #scale_x_continuous(expand=c(0,0.05)) +
            scale_y_continuous(expand=c(0,0.01)) +
            scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
            scale_fill_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
            xlab('Sample size (n)') +
            ylab('Power') +
            ggtitle('MAF=0.02, beta_b=-0.2, beta_g=0.2') +
            theme_classic() +
            theme(legend.position=c(0.8,0.3))
ggsave('figures/power_simulation/n_power.png',n_power,width=4,height=4)

power_plot <- (effect_power + theme(legend.position = 'none')) + (MAF_power + theme(legend.position = 'none')) + n_power

ggsave('figures/power_simulation/power.png',power_plot,width=12,height=4)

