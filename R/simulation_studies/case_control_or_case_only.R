library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
sim_gt_cases <- function(q1,OR1,n_cases){
  f1 <- 1-1/(1+OR1*q1/(1-q1))
  x <- rbinom(n_cases,size=1,prob=f1)
  return(x)
}
compare_cc_co <- function(n_cases,n_controls,q1,OR1_summer,OR1_winter,n_sim){
  y_cc <- c(rep(1,n_cases),rep(0,n_controls))
  y_co <- c(rep(1,n_cases),rep(0,n_cases))
  sim_dat <- lapply(1:n_sim,function(i){
    x_controls <- rbinom(n_controls,size=1,prob=q1)
    x_cases_summer <- sim_gt_cases(q1,OR1_summer,n_cases)
    x_cases_winter <- sim_gt_cases(q1,OR1_winter,n_cases)
    #case control approach
    x_cc_summer <- c(x_cases_summer,x_controls)
    x_cc_winter <- c(x_cases_winter,x_controls)
    mod_cc_summer <- glm(y_cc~x_cc_summer,family='binomial')
    mod_cc_winter <- glm(y_cc~x_cc_winter,family='binomial')
    beta_cc <- mod_cc_summer$coefficients['x_cc_summer'] - mod_cc_winter$coefficients['x_cc_winter']
    se_cc_summer <- summary(mod_cc_summer)$coefficients['x_cc_summer', 2]
    se_cc_winter <- summary(mod_cc_winter)$coefficients['x_cc_winter', 2]
    r_cc <- (n_cases*n_controls/(n_cases + n_controls)) * n_controls/(n_controls^2)
    se_cc <-  sqrt(se_cc_summer^2 + se_cc_winter^2 - 2*se_cc_summer*se_cc_winter*r_cc)
    CI_cc <- beta_cc + c(-1,1)*1.96*se_cc

    #case_only approach
    x_co <- c(x_cases_summer,x_cases_winter)
    mod_co <- glm(y_co~x_co,family='binomial')
    beta_co <- mod_co$coefficients['x_co']
    se_co <- summary(mod_co)$coefficients['x_co', 2]
    CI_co <- beta_co + c(-1,1)*1.96*se_co
    if(beta_cc>10){
      print(table(x_cc_winter))
      print(table(x_cc_summer))
      print(summary(mod_cc_summer))
      print(summary(mod_cc_winter))
    }

    tibble(sim_nr=i,
           method=c('case_control','case_only'),
           beta_est=c(beta_cc,beta_co),
           se=c(se_cc,se_co),
           lower=c(CI_cc[1],CI_co[1]),
           upper=c(CI_cc[2],CI_co[2]))
  }) %>% bind_rows()
  return(sim_dat)
}

visualize_simulation_results <- function(sim_dat,var_name){
  p1 <- group_by(sim_dat,var,method) %>%
        summarise(beta_median=median(beta_est),
                  se_median=median(se)) %>%
        pivot_longer(cols=c('beta_median','se_median'),names_to = 'parameter') %>%
        ggplot(aes(x=factor(var),y=value,fill=method)) +
        geom_col(alpha=0.7,col='black',position='dodge') +
        facet_wrap(~parameter,scales='free') +
        xlab(var_name)

  p2 <- pivot_wider(sim_dat,id_cols=c('var','sim_nr'),names_from='method',values_from=c('beta_est','se')) %>%
        mutate(beta_pct_diff=100*(beta_est_case_control - beta_est_case_only)/beta_est_case_control,
               se_pct_diff=100*(se_case_control-se_case_only)/se_case_control) %>%
        group_by(var) %>%
        summarise(beta_pct_diff.median=median(beta_pct_diff),
                  beta_pct_diff.lower=quantile(beta_pct_diff,0.025),
                  beta_pct_diff.upper=quantile(beta_pct_diff,0.975),
                  se_pct_diff.median=median(se_pct_diff),
                  se_pct_diff.lower=quantile(se_pct_diff,0.025),
                  se_pct_diff.upper=quantile(se_pct_diff,0.975)) %>%
        pivot_longer(cols=-var) %>%
        separate(col=name,into=c('name','metric'),sep='\\.') %>%
        pivot_wider(id_cols=c('var','name'),names_from='metric',values_from='value') %>%
        ggplot(aes(x=as.factor(var))) +
        geom_pointrange(aes(y=median,ymin=lower,ymax=upper)) +
        facet_wrap(~name,scales='free') +
        xlab(var_name) +
        ylab('Percentage difference')

  p1 / p2
}

#vary accross number of cases
n_cases_vec <- c(400,600,800,1000,2000,5000)
case_varying_simulation <- lapply(1:length(n_cases_vec),function(i){
  compare_cc_co(n_cases_vec[i],n_controls=10000,q1=0.1,OR1_summer=1.5,OR1_winter=1.1,n_sim=100) %>%
  mutate(var=n_cases_vec[i])
}) %>% bind_rows()
case_varying_simulation_plot <- visualize_simulation_results(case_varying_simulation,var_name='Number of cases') +
                                plot_annotation(title = 'log(OR_seasonal)=0.31, var_freq=10%')


#vary accross variant frequency
q1_vec <- seq(0.05,0.5,by=0.05)
freq_varying_simulation <- lapply(1:length(q_vec),function(i){
  compare_cc_co(n_cases=1000,n_controls=10000,q1=q1_vec[i],OR1_summer=1.5,OR1_winter=1.1,n_sim=100) %>%
  mutate(var=q1_vec[i])
}) %>% bind_rows()
freq_varying_simulation_plot <- visualize_simulation_results(freq_varying_simulation,var_name='Variant frequency') +
                                plot_annotation(title = 'n_cases=1000, log(OR_seasonal)=0.31')

#Vary accross different seasonal effects
OR1_summer_vec <- seq(1,2,by=0.1)
OR1_summer_varying_simulation <- lapply(1:length(OR1_summer_vec),function(i){
  compare_cc_co(n_cases=1000,n_controls=10000,q1=0.1,OR1_summer=OR1_summer_vec[i],OR1_winter=1.1,n_sim=100) %>%
  mutate(var=round(exp(log(OR1_summer_vec[i])-log(1.1)),2))
}) %>% bind_rows()
OR_summer_varying_simulation_plot <- visualize_simulation_results(OR1_summer_varying_simulation,var_name='Seasonal OR') +
                                      plot_annotation(title = 'n_cases=1000, var_freq=10%')


ggsave('figures/cc_vs_co_case_varying.png',case_varying_simulation_plot)
ggsave('figures/cc_vs_co_freq_varying.png',freq_varying_simulation_plot)
ggsave('figures/cc_vs_co_OR_varying.png',OR_summer_varying_simulation_plot)

