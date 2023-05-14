library(dplyr)
library(tidyr)
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

