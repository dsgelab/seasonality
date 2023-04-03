library(broom)
### DGP of disease endpoint seasonality

dseasonal <- function(t,x_g,beta,phi,a,b){
  ifelse(t>=a & t<b,(1 + (beta[1]+beta[2]*x_g)*sin((2*pi/(b-a)) *(t-a) + phi)) / (b-a),0)
}

Fseasonal <- function(t,x_g,beta,phi,a,b){
  (t-((b-a)/(2*pi))*(beta[1]+beta[2]*x_g)*cos((2*pi/(b-a)) *(t-a) + phi))/(b-a)
}

pseasonal <- function(t,x_g,beta,phi,a,b){
  Fseasonal(t,x_g,beta,phi,a,b) - Fseasonal(a,x_g,beta,phi,a,b)
}

Sseasonal <- function(t,x_g,beta,phi,a,b){
  1-pseasonal(t,x_g,beta,phi,a,b)
}

logit <- function(x){
  log(x/(1-x))
}

inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

inv_beta_map <- function(beta_star,A_inv,b){
  A <- solve(A_inv)
  logit(A %*% (beta_star+b))
}

beta_map <- function(beta,A_inv,b){
  A_inv %*% inv_logit(beta) - b
}

phi_map <- function(phi){
  inv_logit(phi)*pi
}

loglik_seasonal <- function(theta,A_inv,b_t,t,x_g,a,b){
  beta <- beta_map(theta[1:2],A_inv=A_inv,b=b_t)
  phi <- phi_map(theta[3])
  loglik=sum(log(dseasonal(t,x_g,beta,phi,a,b)))
  return(loglik)
}

seasonal_mod <- function(t,x_g,mod_type='additive',a,b){
  if(mod_type=='additive'){
    gamma <- 2
  }else if(mod_type %in% c('recessive','dominant')){
    gamma <- 1
  }else{
    stop('mod_type not recognized')
  }
  A_inv <- matrix(c(2,0,-2/gamma,2/gamma),nrow=2,byrow = T)
  b_t <- c(1,0)
  neg_loglik_seasonal <- function(theta){
    -loglik_seasonal(theta=theta,A_inv=A_inv,b_t=b_t,t=t,x_g=x_g,a=a,b=b)
  }
  opt_res <- optim(par=c(0,0,0),method='BFGS',neg_loglik_seasonal,hessian=T)
  beta_hat <- beta_map(opt_res$par[1:2],A_inv=A_inv,b=b_t)
  phi_hat <- phi_map(opt_res$par[3])
  se_star <- sqrt(diag(solve(opt_res$hessian)))
  conf_interval_star <- cbind(opt_res$par - 1.96*se_star,opt_res$par + 1.96*se_star)
  beta_conf <- beta_map(conf_interval_star[1:2,1:2],A_inv=A_inv,b=b_t)
  phi_conf <- phi_map(conf_interval_star[3,1:2])

  null_beta <- c(inv_beta_map(c(0,beta_hat[2]),A_inv,b_t)[1],
                 inv_beta_map(c(beta_hat[1],0),A_inv,b_t)[2])

  mod_summary <- tibble(param=c('beta_b','beta_g','phi'),
                        est=c(beta_hat,phi_hat),
                        lower=c(beta_conf[,1],phi_conf[1]),
                        upper=c(beta_conf[,2],phi_conf[2]),
                        wald=c(((opt_res$par[1:2]-null_beta)/se_star[1:2])^2,
                               (opt_res$par[3]/se_star[3])^2),
                        pval=pchisq(wald,df=1,lower.tail = F))
  return(mod_summary)
}

get_dist_grid <- function(beta,phi,mod_type='additive',a,b){
  eps <- 1e-6
  t <- seq(a,b-eps,by=0.01)
  gt_mod <- 0:2
  gt_lab <- 0:2
  if(mod_type=='recessive'){
    gt_lab <- gt_mod[1:2]
    gt_mod <- gt_mod[1:2]
  }else if(mod_type=='dominant'){
    gt_lab <- gt_mod[c(1,3)]
    gt_mod <- gt_mod[1:2]
  }
  surv_grid_gt <- lapply(1:length(gt_mod),function(i){
    dens <- dseasonal(t,x_g=gt_mod[i],beta=beta,phi=phi,a=a,b=b)
    cum_dens <- pseasonal(t,x_g=gt_mod[i],beta=beta,phi=phi,a=a,b=b)
    surv <- Sseasonal(t,x_g=gt_mod[i],beta=beta,phi=phi,a=a,b=b)
    tibble(t=t,gt_mod=gt_mod[i],gt_lab=gt_lab[i],dens=dens,cum_dens=cum_dens,surv=surv)
  })
  names(surv_grid_gt) <- gt_lab
  return(surv_grid_gt)
}


simulate_seasonal_dist <- function(n,MAF,beta,phi,mod_type='additive',a,b){
  dist_grid_gt <- get_dist_grid(beta=beta,phi=phi,mod_type=mod_type,a=a,b=b)
  gt_mod <- rbinom(n,size=2,prob=MAF)
  gt_lab <- gt_mod
  if(mod_type=='recessive'){
    gt_lab <- as.integer(gt_mod>0)
    gt_mod <- gt_lab
  }else if(mod_type=='dominant'){
    gt_lab <- 2*as.integer(gt_mod>1)
    gt_mod <- as.integer(gt_mod>1)
  }
  u <- runif(n)
  event_month_dec <- vector('numeric',n)
  for(i in 1:n){
    gt_name <- as.character(gt_lab[i])
    event_month_dec[i] <- dist_grid_gt[[gt_name]]$t[which.min(abs(dist_grid_gt[[gt_name]]$surv-u[i]))]
  }
  return(tibble(gt_mod=gt_mod,gt_lab=gt_lab,EVENT_MONTH_DEC=event_month_dec,surv=u))
}

perturb_surv_times <- function(surv_times,mu,sigma,a,b){
  n <- length(surv_times)
  diagnosis_lag <- rgamma(n,shape=(mu/sigma)^2,rate=mu/sigma^2)
  anonymization_perturb <- sample(seq_len(15),size = n,replace=T)/28 #convert to months
  anonymization_perturb <- ifelse(runif(n)<0.5,-1*anonymization_perturb,anonymization_perturb)
  surv_times_perturb <- (surv_times+diagnosis_lag+anonymization_perturb - a) %% (b-a) + a
  return(surv_times_perturb)
}

run_seasonality_gam_sim <- function(dat,a,b){
  k_seasonal <- 6
  f_seasonal <- 'count ~ s(EVENT_MONTH, k=k_seasonal, bs="cp")'
  mod_seasonal <- gam(as.formula(f_seasonal), family=quasipoisson(), data=dat, knots=list(EVENT_MONTH = c(a-0.5, b-0.5)), scale=-1)
  return(mod_seasonal)
}

get_seasonal_phenotype <- function(pheno_dat,endpoint_id,seasonal_spline,a,b){
  pheno_dat[,endpoint_id,drop=T] <- ifelse(pheno_dat[,endpoint_id,drop=T]>(b-0.5),pheno_dat[,endpoint_id,drop=T]-(b-a),pheno_dat[,endpoint_id,drop=T])
  seasonal_pheno_dat <- mutate(pheno_dat,
                               seasonal_val=approx(x=seasonal_spline$month,
                                                   y=seasonal_spline$seasonal_val,
                                                   xout=pheno_dat[,endpoint_id,drop=T])$y,
                               seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                               seasonal_val_binary=as.integer(seasonal_val > 0))
  return(seasonal_pheno_dat)
}

run_seasonal_models <- function(monthly_counts,pheno_dat,endpoint_id,a,b){
  mod_dat <- filter(monthly_counts,ENDPOINT==endpoint_id)
  mod_counts <- run_seasonality_gam_sim(mod_dat,a=a,b=b)
  months_vec <- seq(a-0.5,b-0.5,by=0.01)
  seasonal_pred <- predict(mod_counts,newdata=data.frame(EVENT_MONTH=months_vec),type='terms')[,'s(EVENT_MONTH)']
  seasonal_spline <- tibble(month=months_vec,seasonal_val=seasonal_pred)
  seasonal_pheno <- get_seasonal_phenotype(pheno_dat=pheno_dat,endpoint_id=endpoint_id,seasonal_spline = seasonal_spline,a=a,b=b) %>%
                    mutate(seasonal_val_01=(seasonal_val-min(seasonal_val))/(max(seasonal_val)-min(seasonal_val)))
  binary_mod <- glm(seasonal_val_binary ~ gt_mod,data=seasonal_pheno,family='binomial')
  qt_raw_mod <- lm(seasonal_val_01 ~ gt_mod,data=seasonal_pheno)
  qt_mod <- lm(seasonal_val_qt ~ gt_mod,data=seasonal_pheno)
  return(list(seasonal_pheno=seasonal_pheno,binary_mod=binary_mod,qt_raw_mod=qt_raw_mod,qt_mod=qt_mod))
}


sim_pipeline <- function(n,MAF,beta,phi,mod_type='additive',a,b,mu_lag,sigma_lag){
  sim_dat <- simulate_seasonal_dist(n=n,MAF=MAF,beta=beta,phi=phi,mod_type=mod_type,a=a,b=b)
  sim_dat$EVENT_MONTH_DEC_PB <- perturb_surv_times(surv_times=sim_dat$EVENT_MONTH_DEC,mu=mu_lag,sigma=sigma_lag,a=a,b=b)
  monthly_counts <- mutate(sim_dat,month=floor(EVENT_MONTH_DEC),
                           month_pb=floor(EVENT_MONTH_DEC_PB)) %>%
                    pivot_longer(cols=c('month','month_pb'),names_to = 'ENDPOINT',values_to = 'EVENT_MONTH') %>%
                    group_by(ENDPOINT,EVENT_MONTH) %>%
                    summarise(count=length(EVENT_MONTH),.groups='drop') %>%
                    mutate(ENDPOINT=as.character(factor(ENDPOINT,levels=c('month','month_pb'),
                                                        labels=c('EVENT_MONTH_DEC','EVENT_MONTH_DEC_PB'))),
                           ENDPOINT_LAB=as.character(factor(ENDPOINT,levels=c('EVENT_MONTH_DEC','EVENT_MONTH_DEC_PB'),
                                                        labels=c('True disease onset','Perturbed disease diagnosis'))))
  seasonal_mods <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC',a=a,b=b)
  seasonal_mods_pb <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC_PB',a=a,b=b)
  mod <- c('binary_mod','qt_raw_mod','qt_mod')
  mod_summary <- lapply(mod,function(m){
    bind_rows(tidy(seasonal_mods[[m]]) %>% mutate(perturbed=F),
              tidy(seasonal_mods_pb[[m]]) %>% mutate(perturbed=T)) %>%
    mutate(mod=m) %>% filter(term=='gt_mod') %>% select(-term)
  }) %>% bind_rows()
  return(list(sim_dat=sim_dat,mod_summary=mod_summary))
}

plot_surv_grid_panel <- function(surv_grid){
  ggplot(surv_grid,aes(t,value,col=as.factor(gt))) +
    geom_line(size=1) +
    facet_wrap(sim~name,ncol=2,scales='free') +
    scale_x_continuous(breaks=seq_len(12)) +
    scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='Genotype') +
    xlab('Month') +
    ylab('Log hazard (left panel) and Survival (right panel)') +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_blank())
}

get_seasonal_comp <- function(dat,mod_gam,a,b){
  months_vec <- seq(a-0.5,b-0.5,by=0.1)
  seasonal_pred <- predict(mod_gam,newdata=get_grid(dat,type='seasonal',seasonal_spline_avg=NULL),type='terms',se.fit=T)
  seasonal_pred_dat <- tibble(month=months_vec,
                              est=seasonal_pred$fit[,'s(EVENT_MONTH)'],
                              lower=est - 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'],
                              upper=est + 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'])
  return(seasonal_pred_dat)
}



## Plots
seasonality_plot_sim <- function(seasonal_pred_dat,a,b){
  ggplot(seasonal_pred_dat) +
    geom_line(aes(x=month,
                  y=est,
                  col=ENDPOINT_LAB)) +
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=ENDPOINT_LAB),alpha=0.3) +
    geom_hline(yintercept=0,linetype='dashed',color='red') +
    scale_x_continuous(breaks=seq(a,b)) +
    scale_fill_manual(values=c("#0072B5FF","#E18727FF"),name='') +
    scale_color_manual(values=c("#0072B5FF","#E18727FF"),name='') +
    xlab('Month number') +
    ylab('Smooth term') +
    ggtitle('Seasonal component') +
    theme_bw() +
    theme(axis.title.x=element_text(size=14),
          axis.text.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          axis.text.y=element_text(size=14))
}


KM_gg <- function(survfit_obj,title=''){
  cols <- c("#BC3C29FF","#0072B5FF","#E18727FF")
  gt_strata <- as.integer(paste0(gsub('gt_lab=','',names(survfit_obj$strata))))
  cols <- cols[gt_strata+1]
  survfit_dat <- tibble(time=rep(survfit_obj$time,survfit_obj$n.event),
                        surv=rep(survfit_obj$surv,survfit_obj$n.event),
                        strata=rep(paste0(gt_strata,' (n=',survfit_obj$n,')'),survfit_obj$n)) %>%
                  distinct()
  ggplot(survfit_dat,aes(time,surv,col=factor(strata))) +
    geom_step() +
    scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
    scale_color_manual(values=cols,name='Genotype') +
    ggtitle(title) +
    xlab('Time in months') +
    ylab('Survival') +
    theme_classic() +
    theme(legend.position=c(0.8,0.8))
}

### DGP for
