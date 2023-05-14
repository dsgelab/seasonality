library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(survival)
library(mgcv)
library(patchwork)
library(data.table)
library(ggthemes)
source('seasonality_dashboard/utils/process_underlying_seasonality_utils.R')
source('seasonality_dashboard/utils/seasonality_utils.R')
source('aux/tix_plotting.R')

###### Theoretical genetical seasonal effect ######
# Functions to produce figures demonstrating amplitude and phase seasonal effects
dgp_mu <- function(beta,phi,x_g,t,P){
  beta['g']*x_g + (beta['a'] + beta['ag']*x_g)*sin(t*2*pi/P + phi + beta['pg']*x_g)
}

get_mu_grid <- function(beta,phi,mod_type,P){
  t_seq <- seq(0,P,by=0.01)
  gt_mod <- 0:2
  gt_lab <- 0:2
  if(mod_type=='recessive'){
    gt_lab <- gt_mod[1:2]
    gt_mod <- gt_mod[1:2]
  }else if(mod_type=='dominant'){
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

#Theoretical genetic seasonal effect curves
input <- list(P=12, # period
              beta=list('g'=rep(0,6),  # additive genetic effect
                     'a'=rep(0.15,6),  # baseline seasonal/circadian amplitude
                     'ag'=c(0.1,0,0.1,-0.1,0,-0.1),   # genetic-season amplitude interaction
                     'pg'=c(0,0.3,0.3,0,-0.3,-0.3)), # genetic-season phase interaction
              phi=1, # baseline phase
              mod_type='additive')

plot_idx <- c(1,3,5,2,4,6)
plot_dat <- lapply(1:length(input$beta$ag),function(i){
              mu_grid <- get_mu_grid(beta=sapply(input$beta,function(b) b[i]),phi=input$phi,mod_type=input$mod_type,P=input$P)
              mutate(mu_grid,t=t+1,
                             idx=toupper(letters[plot_idx[i]]))
            }) %>% bind_rows()

cols <- c("#BC3C29FF","#0072B5FF","#E18727FF")
combined_plot <- ggplot(data=plot_dat) +
                geom_line(aes(t,mu,col=factor(gt_lab)),linewidth=0.8) +
                facet_wrap(~factor(idx,levels=toupper(letters[plot_idx])),ncol=3) +
                geom_hline(yintercept=0,linetype='dashed') +
                scale_color_manual(values=cols,name='Genotype') +
                scale_x_continuous(breaks=seq(1,input$P)) +
                ylab('') +
                xlab('Month') +
                theme_classic() +
                annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,linewidth=0.8)+
                annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,linewidth=0.8) +
                theme(strip.text = element_text(size = 12, margin = margin(),hjust = 0),
                      strip.background = element_blank(),
                      axis.text.x=element_text(size=7))
                      #panel.border = element_rect(colour = "black", fill = NA))

save_tikz_plot(combined_plot,width=7,height=4,filename = 'tex/genetic_seasonal_effect_theory.tex')

###### Statistical inference of seasonality density - Demonstration ######
# Mapping from parallelogram to plane for statistical inference of DGP model
gamma <- 2
A_inv <- matrix(c(2,0,-2/gamma,2/gamma),nrow=2,byrow = T)
b_t <- c(1,0)
beta <- c(-0.5,0.5)
beta_square <- solve(A_inv) %*% (beta+b_t)
beta_plane <- logit(beta_square)
point_dat <- tibble()
parallelogram_plot <- tibble(x=c(-1,-1,1,1),y=c(0,1,0,-1),xend=c(-1,1,1,-1),yend=c(1,0,-1,0)) %>%
                      ggplot() +
                      geom_point(x=beta[1],y=beta[2],col='red',size=2) +
                      geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
                      geom_vline(xintercept=0,linetype='dashed') +
                      geom_hline(yintercept=0,linetype='dashed') +
                      scale_x_continuous(limits=c(-2,2)) +
                      scale_y_continuous(limits=c(-2,2)) +
                      xlab('$beta_b$') +
                      ylab('$beta_{g}$') +
                      theme_classic()

rectangle_plot <- tibble(x=c(0,0,1,1),y=c(0,1,1,0),xend=c(0,1,1,0),yend=c(1,1,0,0)) %>%
                  ggplot() +
                  geom_point(x=beta_square[1],y=beta_square[2],col='red',size=2) +
                  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
                  geom_vline(xintercept=0,linetype='dashed') +
                  geom_hline(yintercept=0,linetype='dashed') +
                  scale_x_continuous(limits=c(-2,2)) +
                  scale_y_continuous(limits=c(-2,2)) +
                  xlab(expression('$beta_b$')) +
                  ylab('') +
                  theme_classic()
plane_plot <- tibble(x=c(0,0,1,1),y=c(0,1,1,0),xend=c(0,1,1,0),yend=c(1,1,0,0)) %>%
              ggplot() +
              geom_point(x=beta_plane[1],y=beta_plane[2],col='red',size=2) +
              geom_vline(xintercept=0,linetype='dashed') +
              geom_hline(yintercept=0,linetype='dashed') +
              scale_x_continuous(limits=c(-2,2)) +
              scale_y_continuous(limits=c(-2,2)) +
              xlab('$beta_b$') +
              ylab('') +
              theme_classic()
combined_mapping_plot <- parallelogram_plot + rectangle_plot + plane_plot
save_tikz_plot(combined_mapping_plot,width=7,height=7/3,filename = 'tex/DGP_model_mapping.tex')

###### Simulation based approach to validate seasonality GWAS framework ######
a <- 1;b <- 13
models <- c('binary_mod','qt_raw_mod')
model_names <- c('Binary phenotype','Quantitative phenotype')
cols <- c("#BC3C29FF","#0072B5FF","#E18727FF")
cols_mod <- cols[1:length(models)]
beta_int_seq <- seq(-0.4,0.4,by=0.05)

###### Correlation between true effect and seasonality GWAS effects ######
#R2 estimates from effect sizes
genetic_effect_simulation <- read_tsv('results/simulations/genetic_effect_simulation.tsv')
lapply(1:max(genetic_effect_simulation$sim_nr), function(i) {
  mod_dat <- filter(genetic_effect_simulation,sim_nr==i)
  mod_table_binary <- lm(beta_int ~ estimate,data=filter(mod_dat,mod=='binary_mod')) %>%
    glance(mod) %>% mutate(sim_nr=i,mod='binary') %>% select(sim_nr,mod,everything())
  mod_table_qt <- lm(beta_int ~ estimate,data=filter(mod_dat,mod=='qt_raw_mod')) %>%
    glance(mod) %>% mutate(sim_nr=i,mod='qt') %>% select(sim_nr,mod,everything())
  bind_rows(mod_table_binary,mod_table_qt)
}) %>% bind_rows() %>% group_by(mod) %>% summarise(mean_r2=mean(r.squared),lower=quantile(r.squared,probs = 0.025),upper=quantile(r.squared,probs = 0.975))

###### Demonstration of simulation procedure ######
#Simulate DGP and run simulation procedure
set.seed(100)
beta_b <- 0.3;beta_ag <- -0.2;phi <- pi/2
sim_dat <- simulate_seasonal_dist(n=40000,MAF=0.1,beta=c(beta_b,beta_ag),phi=phi,mod_type='additive',a=a,b=b)
dist_grid <- get_dist_grid(beta=c(beta_b,beta_ag),phi=phi,mod_type='additive',a=a,b=b) %>% bind_rows()
dens_plot <- ggplot(dist_grid,aes(t,dens,col=as.factor(gt_lab))) +
  geom_line(size=1) +
  scale_x_continuous(breaks=seq_len(12),expand=c(0,0.2)) +
  scale_y_continuous(limits=c(0,2/12)) +
  scale_color_manual(values=cols,name='Genotype') +
  xlab('Month number ') +
  ylab('Density') +
  theme_classic() +
  #theme_custom +
  theme(legend.position='none')

surv_obj_sim <- Surv(time=sim_dat$EVENT_MONTH_DEC,event=rep(1,nrow(sim_dat)))
survfit_sim <- survfit(surv_obj_sim ~ gt_lab,data=sim_dat)
KM_plot <- KM_gg(survfit_sim,'') + scale_x_continuous(breaks=seq_len(12)) + xlab('Month number')
monthly_counts <- mutate(sim_dat,EVENT_MONTH=floor(EVENT_MONTH_DEC)) %>%
  mutate(ENDPOINT='EVENT_MONTH_DEC',
         ENDPOINT_LAB='True disease onset') %>%
  group_by(ENDPOINT,ENDPOINT_LAB,EVENT_MONTH) %>%
  summarise(count=length(EVENT_MONTH))
monthly_counts_plot <- ggplot(monthly_counts,aes(EVENT_MONTH,count)) +
  geom_col(position='dodge',fill="#0072B5FF") +
  scale_x_continuous(breaks=seq(1,12,by=2),expand=c(0,0.2)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('Month number') +
  ylab('Count') +
  theme_classic()

mod_counts_true <- run_seasonality_gam_sim(monthly_counts,a=a,b=b)

seasonal_comp <- get_seasonal_comp(mod_dat_true,mod_counts_true,a=a,b=b) %>% mutate(ENDPOINT_LAB='True disease onset')
GAM_seasonality_plot <- seasonality_plot_sim(seasonal_comp,a=a,b=b) + scale_x_continuous(breaks=seq(1,12,by=2)) + theme_classic() + theme(legend.position = 'none') + ggtitle('')#  + theme_custom

months_vec <- seq(a-0.5,b-0.5,by=0.01)
seasonal_pred <- predict(mod_counts_true,newdata=data.frame(EVENT_MONTH=months_vec),type='terms')[,'s(EVENT_MONTH)']
seasonal_spline <- tibble(month=months_vec,seasonal_val=seasonal_pred)
seasonal_pheno <- get_seasonal_phenotype(pheno_dat=sim_dat,endpoint_id='EVENT_MONTH_DEC',seasonal_spline = seasonal_spline,a=a,b=b) %>%
  mutate(seasonal_val_01=(seasonal_val-min(seasonal_val))/(max(seasonal_val)-min(seasonal_val)))

phenotype_plot <- ggplot(seasonal_pheno) +
  geom_histogram(aes(seasonal_val_01,fill=factor(seasonal_val_binary)),bins=80) +
  scale_fill_manual(values=c("#E18727FF","#0072B5FF"),name='Binary phenotype') +
  scale_x_continuous(expand=c(0.005,0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab('Continuous phenotype') +
  ylab('Count') +
  theme_classic() +
  theme(legend.position=c(0.45,0.5))


sim_plot <- ((dens_plot + ggtitle('A')) + (KM_plot + ggtitle('B'))) / ((monthly_counts_plot + ggtitle('C')) + (GAM_seasonality_plot + ggtitle('D')) + (phenotype_plot + ggtitle('E')))
save_tikz_plot(sim_plot,width=7,height=6,filename = 'tex/simulation_overview_plot.tex')


###### Power plot + diagnosis lag sensitivity ######
mu_lag_seq <- seq(0,6,by=0.5)
lag_simulation <- read_tsv('results/simulations/lag_simulation.tsv')

pval_lag_plot <- filter(lag_simulation,perturbed) %>%
                filter(mod %in% models) %>%
                mutate(mod=factor(mod,levels=models,labels=model_names)) %>%
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
                  geom_hline(yintercept=-log10(5e-8),linetype='dashed') +
                  #ggtitle('P-value of associations on perturbed disease onset data') +
                  scale_color_manual(values=cols_mod,name='') +
                  xlab('Mean diagnosis lag (in months)') +
                  ylab(expression(-log[10](p))) +
                  theme_bw()


## Power simulation
# vary effect size
n_sim <- 100
effect_power_simulation <- read_tsv('results/simulations/effect_power_simulation')

effect_power <- filter(effect_power_simulation,!perturbed) %>%
                filter(mod %in% models) %>%
                mutate(mod=factor(mod,levels=models,labels=model_names)) %>% group_by(mod,beta_int) %>%
                summarise(num_GWAS_signif=sum(p.value<5e-8),
                          power=num_GWAS_signif/n_sim,
                          power_se=sqrt(power*(1-power)/n_sim)) %>%
                ggplot() +
                geom_line(aes(beta_int,power,col=mod)) +
                geom_ribbon(aes(beta_int,ymin=power-1.96*power_se,ymax=power+1.96*power_se,fill=mod),alpha=0.1) +
                scale_x_continuous(expand=c(0,0.01)) +
                scale_y_continuous(expand=c(0,0.01)) +
                scale_color_manual(values=cols_mod,name='') +
                scale_fill_manual(values=cols_mod,name='') +
                xlab('$beta_g$') +
                ylab('Power') +
                ggtitle('A MAF=2\\%, $beta_b$=-0.2, n=40000') +
                theme_classic() +
                theme(legend.position=c(0.7,0.3))

#vary MAF
MAF_power_simulation <- read_tsv('results/simulations/MAF_power_simulation.tsv')

MAF_power <- filter(MAF_power_simulation,!perturbed) %>%
            filter(mod %in% models) %>%
            mutate(mod=factor(mod,levels=models,labels=model_names),
                   maf=maf*100) %>%
            group_by(mod,maf) %>%
            summarise(num_GWAS_signif=sum(p.value<5e-8),
                      power=num_GWAS_signif/n_sim,
                      power_se=sqrt(power*(1-power)/n_sim)) %>%
            ggplot() +
            geom_line(aes(maf,power,col=mod)) +
            geom_ribbon(aes(maf,ymin=power-1.96*power_se,ymax=power+1.96*power_se,fill=mod),alpha=0.1) +
            scale_x_continuous(breaks=seq(0,30,by=5),expand=c(0,0.01)) +
            scale_y_continuous(expand=c(0,0.01)) +
            scale_color_manual(values=cols_mod,name='') +
            scale_fill_manual(values=cols_mod,name='') +
            xlab('MAF (\\%)') +
            ylab('Power') +
            ggtitle('B $beta_b$=-0.4, $beta_g$=0.1, n=40000') +
            theme_classic() +
            theme(legend.position=c(0.7,0.3))

#vary sample size
n_power_simulation <- read_tsv('results/simulations/n_power_simulation.tsv')

n_power <- filter(n_power_simulation,!perturbed) %>%
            filter(mod %in% models) %>%
            mutate(mod=factor(mod,levels=models,labels=model_names)) %>%
            group_by(mod,n) %>%
            summarise(num_GWAS_signif=sum(p.value<5e-8),
                      power=num_GWAS_signif/n_sim,
                      power_se=sqrt(power*(1-power)/n_sim)) %>%
            ggplot() +
            geom_line(aes(n,power,col=mod)) +
            geom_ribbon(aes(n,ymin=power-1.96*power_se,ymax=power+1.96*power_se,fill=mod),alpha=0.1) +
            scale_y_continuous(expand=c(0,0.01)) +
            scale_color_manual(values=cols_mod,name='') +
            scale_fill_manual(values=cols_mod,name='') +
            xlab('Sample size (n)') +
            ylab('Power') +
            ggtitle('C MAF=2\\%, $beta_b$=-0.2, $beta_g$=0.2') +
            theme_classic() +
            theme(legend.position=c(0.7,0.3))

power_plot <- ((effect_power + theme(legend.position = 'none')) + (MAF_power)) /
              ((n_power + theme(legend.position = 'none')) + (pval_lag_plot + ggtitle('D') + theme(legend.position = 'none')))
save_tikz_plot(power_plot,width=7.5,height=6,filename = 'tex/simulation_power.tex')


# baseline_effect_simulation <- read_tsv('results/simulations/baseline_effect_simulation.tsv')
#
# effect_true_onset_baseline_plot <- filter(baseline_effect_simulation,!perturbed) %>%
#   group_by(mod,beta_base) %>%
#   summarise(avg_effect=mean(estimate),
#             lower_effect=quantile(estimate,0.025),
#             upper_effect=quantile(estimate,0.975)) %>%
#   ggplot(aes(x=beta_base,
#              y=avg_effect,
#              ymin=lower_effect,
#              ymax=upper_effect,
#              col=factor(mod))) +
#   geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
#   scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
#   xlab(expression(beta[b])) +
#   ylab(expression(hat(beta)[pheno])) +
#   theme_bw()
#
#
# pval_true_onset_baseline_plot <- filter(baseline_effect_simulation,!perturbed) %>%
#   group_by(mod,beta_base) %>%
#   mutate(minus_log10_p=-log10(p.value)) %>%
#   summarise(avg_effect=mean(minus_log10_p),
#             lower_effect=quantile(minus_log10_p,0.025),
#             upper_effect=quantile(minus_log10_p,0.975)) %>%
#   ggplot(aes(x=beta_base,
#              y=avg_effect,
#              ymin=lower_effect,
#              ymax=upper_effect,
#              col=factor(mod))) +
#   geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
#   geom_abline(intercept=-log10(5e-8),slope=0,linetype='dashed') +
#   scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='') +
#   xlab(expression(beta[b])) +
#   ylab(expression(-log[10](p))) +
#   theme_bw()
#
# true_onset_baseline_plot <- (effect_true_onset_baseline_plot + theme(legend.position='none')) + pval_true_onset_baseline_plot
# ggsave('figures/power_simulation/effect_pval_true_onset_baseline.png',true_onset_baseline_plot,width=8,height=4)

# effect_true_onset_plot <- filter(genetic_effect_simulation,!perturbed) %>%
#   filter(mod %in% models) %>%
#   mutate(mod=factor(mod,levels=models,labels=model_names)) %>%
#   group_by(mod,beta_int) %>%
#   summarise(avg_effect=mean(estimate),
#             lower_effect=quantile(estimate,0.025),
#             upper_effect=quantile(estimate,0.975)) %>%
#   ggplot(aes(x=beta_int,
#              y=avg_effect,
#              ymin=lower_effect,
#              ymax=upper_effect,
#              col=factor(mod))) +
#   geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
#   geom_abline(intercept=0,slope=1) +
#   scale_color_manual(values=cols_mod,name='') +
#   xlab(expression(beta[g])) +
#   ylab(expression(hat(beta)[pheno])) +
#   theme_bw()
#
#
# pval_true_onset_plot <- filter(genetic_effect_simulation,!perturbed) %>%
#   filter(mod %in% models) %>%
#   mutate(mod=factor(mod,levels=models,labels=model_names)) %>%
#   group_by(mod,beta_int) %>%
#   mutate(minus_log10_p=-log10(p.value)) %>%
#   summarise(avg_effect=mean(minus_log10_p),
#             lower_effect=quantile(minus_log10_p,0.025),
#             upper_effect=quantile(minus_log10_p,0.975)) %>%
#   ggplot(aes(x=beta_int,
#              y=avg_effect,
#              ymin=lower_effect,
#              ymax=upper_effect,
#              col=factor(mod))) +
#   geom_pointrange(position = position_jitterdodge(dodge.width = 0.01, jitter.width = 0)) +
#   geom_hline(yintercept=-log10(5e-8),linetype='dashed') +
#   scale_color_manual(values=cols_mod,name='') +
#   xlab(expression(beta[g])) +
#   ylab(expression(-log[10](p))) +
#   theme_bw()
#
# true_onset_plot <- (effect_true_onset_plot + theme(legend.position='none')) + pval_true_onset_plot


# effect_lag_plot <- filter(lag_simulation,perturbed) %>%
#   filter(mod %in% models) %>%
#   mutate(mod=factor(mod,levels=models,labels=model_names)) %>%
#   group_by(mod,mu_lag) %>%
#   summarise(avg_effect=mean(estimate),
#             lower_effect=quantile(estimate,0.025),
#             upper_effect=quantile(estimate,0.975)) %>%
#   ggplot(aes(x=mu_lag,
#              y=avg_effect,
#              ymin=lower_effect,
#              ymax=upper_effect,
#              col=factor(mod))) +
#   geom_pointrange(size=0.7,position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0)) +
#   geom_abline(intercept=0,slope=0,linetype='longdash') +
#   #ggtitle('Effect size for perturbed disease onset data') +
#   #scale_x_continuous(breaks=round(beta_int_seq,2)) +
#   scale_color_manual(values=cols_mod,name='') +
#   xlab('Mean diagnosis lag (in months)') +
#   ylab(expression(hat(beta)[pheno])) +
#   theme_bw()

# seasonal_mods <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC',a=a,b=b)
#
# summary_dat <- bind_rows(list(binary_mod=tidy(seasonal_mods$binary_mod),
#                               qt_mod=tidy(seasonal_mods$qt_mod)),.id='mod_type') %>%
#   filter(term=='gt_mod') %>%
#   mutate(mod_type=gsub('_mod','',mod_type)) %>%
#   separate(mod_type,into=c('model_type','data_type')) %>%
#   mutate(data_type=ifelse(is.na(data_type),'True disease onset','Perturbed disease diagnosis'))
# effect_plot <- ggplot(summary_dat) +
#   geom_pointrange(aes(x=model_type,y=estimate,ymin=estimate-1.96*std.error,ymax=estimate+1.96*std.error,col=factor(data_type,levels=c('True disease onset','Perturbed disease diagnosis'))),
#                   position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0),size=1) +
#   scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
#   xlab('Model') +
#   ylab('Effect size') +
#   theme_classic() +
#   theme_custom
#
# pval_plot <- ggplot(summary_dat) +
#   geom_col(aes(x=model_type,y=-log10(p.value),fill=factor(data_type,levels=c('True disease onset','Perturbed disease diagnosis'))),position = 'dodge',width=0.6) +
#   geom_abline(slope=0,intercept=-log10(5e-8),color='red',linetype='dashed') +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_fill_manual(values=c("#0072B5FF","#E18727FF"),name='') +
#   xlab('Model') +
#   ylab(expression(-log10(pval))) +
#   theme_classic() +
#   theme_custom



