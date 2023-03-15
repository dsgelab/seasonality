# log_hseasonal <- function(t,gt,beta,phi,s,P,baseline=T){
#   log_h <- beta[1] + beta[2]*gt + beta[3]*sin(2*pi*(t-s) / P + phi) + beta[4]*gt*sin(2*pi*(t-s)/P + phi)
#   if(baseline){
#     log_h <- log_h - log(P+s-t)
#   }
#   return(log_h)
# }
#
# hseasonal_opt <- function(t,gt,beta,phi,s,P,baseline=T){
#   h <- beta[1] + beta[2]*gt + beta[3]*sin(2*pi*(t-s) / P + phi) + beta[4]*gt*sin(2*pi*(t-s)/P + phi)
#   if(baseline){
#     h <- h/(P+s-t)
#   }
#   return(h)
# }
#
#
# get_surv_grid <- function(beta,phi,baseline=T,P=12,s=1){
#   eps <- 1e-6
#   t <- seq(s,P+s-eps,by=0.01)
#   gt <- 0:2
#   hseasonal <- function(t,gt){
#     exp(log_hseasonal(t=t,gt=gt,beta=beta,phi=phi,s=s,P=P,baseline=baseline))
#   }
#   surv_grid_gt <- lapply(1:length(gt),function(i){
#     log_haz <- log_hseasonal(t=t,gt=gt[i],beta=beta,phi=phi,s=s,P=P,baseline=baseline)
#     haz <- exp(log_haz)
#     #haz <- hseasonal_opt(t=t,gt=gt[i],beta=beta,phi=phi,s=s,P=P,baseline=baseline)
#     cum_haz <- sapply(t,function(u) integrate(hseasonal,lower=1,upper=u,gt=gt[i])$value)
#     surv <- exp(-cum_haz)
#     cum_dens <- 1-surv
#     dens <- haz*surv
#     tibble(t=t,gt=gt[i],log_haz=log_haz,haz=haz,cum_haz=cum_haz,surv=surv,cum_dens=cum_dens,dens=dens)
#   })
# }
#
# simulate_seasonal_surv <- function(n,MAF,beta,phi,P=12,s=1){
#   surv_grid_gt <- get_surv_grid(beta=beta,phi=phi,P=P,s=s)
#   n <- 40000
#   MAF <- 0.05
#   gt <- rbinom(n,size=2,prob=MAF)
#   u <- runif(n)
#   event_month_dec <- vector('numeric',n)
#   for(i in 1:n){
#     event_month_dec[i] <- surv_grid_gt[[gt[i]+1]]$t[which.min(abs(surv_grid_gt[[gt[i]+1]]$surv-u[i]))]
#   }
#   return(tibble(gt=gt,EVENT_MONTH_DEC=event_month_dec,surv=u))
# }

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

loglik_seasonal <- function(theta){
  loglik=sum(log(dseasonal(t,x_g,theta[1:2],theta[3],a,b)))
  return(loglik)
}

get_dist_grid <- function(beta,phi,a=0,b=12){
  eps <- 1e-6
  t <- seq(a,b-eps,by=0.01)
  gt <- 0:2
  surv_grid_gt <- lapply(1:length(gt),function(i){
    dens <- dseasonal(t,x_g=gt[i],beta=beta,phi=phi,a=a,b=b)
    cum_dens <- pseasonal(t,x_g=gt[i],beta=beta,phi=phi,a=a,b=b)
    surv <- Sseasonal(t,x_g=gt[i],beta=beta,phi=phi,a=a,b=b)
    tibble(t=t,gt=gt[i],dens=dens,cum_dens=cum_dens,surv=surv)
  })
}


simulate_seasonal_dist <- function(n,MAF,beta,phi,a=0,b=12){
  dist_grid_gt <- get_dist_grid(beta=beta,phi=phi,a=a,b=b)
  gt <- rbinom(n,size=2,prob=MAF)
  u <- runif(n)
  event_month_dec <- vector('numeric',n)
  for(i in 1:n){
    event_month_dec[i] <- dist_grid_gt[[gt[i]+1]]$t[which.min(abs(dist_grid_gt[[gt[i]+1]]$surv-u[i]))]
  }
  return(tibble(gt=gt,EVENT_MONTH_DEC=event_month_dec,surv=u))
}

perturb_surv_times <- function(surv_times,mu,sigma,P=12,s=1){
  n <- length(surv_times)
  diagnosis_lag <- rgamma(n,shape=(mu/sigma)^2,rate=mu/sigma^2)
  anonymization_perturb <- sample(seq_len(15),size = n,replace=T)/28 #convert to months
  anonymization_perturb <- ifelse(runif(n)<0.5,-1*anonymization_perturb,anonymization_perturb)
  surv_times_perturb <- (surv_times+diagnosis_lag+anonymization_perturb - s) %% P + s
  return(surv_times_perturb)
}

run_seasonality_gam_sim <- function(dat){
  k_seasonal <- 6
  f_seasonal <- 'count ~ s(EVENT_MONTH, k=k_seasonal, bs="cp")'
  mod_seasonal <- gam(as.formula(f_seasonal), family=nb(), data=dat, knots=list(EVENT_MONTH = c(0.5, 12.5)), scale=-1)
  return(mod_seasonal)
}

get_seasonal_phenotype <- function(pheno_dat,endpoint_id,seasonal_spline){
  seasonal_pheno_dat <- mutate(pheno_dat,
                               seasonal_val=approx(x=seasonal_spline$month,
                                                   y=seasonal_spline$seasonal_val,
                                                   xout=pheno_dat[,endpoint_id,drop=T])$y,
                               seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                               seasonal_val_binary=as.integer(seasonal_val > 0))
  return(seasonal_pheno_dat)
}

run_seasonal_models <- function(monthly_counts,pheno_dat,endpoint_id){
  mod_dat <- filter(monthly_counts,ENDPOINT==endpoint_id)
  mod_counts <- run_seasonality_gam_sim(mod_dat)
  months_vec <- seq(0.5,12.5,by=0.01)
  seasonal_pred <- predict(mod_counts,newdata=data.frame(EVENT_MONTH=months_vec),type='terms')[,'s(EVENT_MONTH)']
  seasonal_spline <- tibble(month=months_vec,seasonal_val=seasonal_pred)
  seasonal_pheno <- get_seasonal_phenotype(pheno_dat=pheno_dat,endpoint_id=endpoint_id,seasonal_spline = seasonal_spline)
  binary_mod <- glm(seasonal_val_binary ~ gt,data=seasonal_pheno,family='binomial')
  qt_mod <- lm(seasonal_val_qt ~ gt,data=seasonal_pheno)
  return(list(seasonal_pheno=seasonal_pheno,binary_mod=binary_mod,qt_mod=qt_mod))
}

sim_pipeline <- function(n,MAF,beta,phi,a,b,mu_lag,sigma_lag){
  sim_dat <- simulate_seasonal_dist(n=n,MAF=MAF,beta=beta,phi=phi,a=a,b=b)
  sim_dat$EVENT_MONTH_DEC_PB <- perturb_surv_times(surv_times=sim_dat$EVENT_MONTH_DEC,mu=mu_lag,sigma=sigma_lag)
  monthly_counts <- mutate(sim_dat,month=floor(EVENT_MONTH_DEC),
                           month_pb=floor(EVENT_MONTH_DEC_PB)) %>%
                    pivot_longer(cols=c('month','month_pb'),names_to = 'ENDPOINT',values_to = 'EVENT_MONTH') %>%
                    group_by(ENDPOINT,EVENT_MONTH) %>%
                    summarise(count=length(EVENT_MONTH),.groups='drop') %>%
                    mutate(ENDPOINT=as.character(factor(ENDPOINT,levels=c('month','month_pb'),
                                                        labels=c('EVENT_MONTH_DEC','EVENT_MONTH_DEC_PB'))),
                           ENDPOINT_LAB=as.character(factor(ENDPOINT,levels=c('EVENT_MONTH_DEC','EVENT_MONTH_DEC_PB'),
                                                        labels=c('True disease onset','Perturbed disease diagnosis'))))
  seasonal_mods <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC')
  seasonal_mods_pb <- run_seasonal_models(monthly_counts,sim_dat,endpoint_id='EVENT_MONTH_DEC_PB')
  return(list(sim_dat=sim_dat,
              binary_mod=seasonal_mods$binary_mod,
              qt_mod=seasonal_mods$qt_mod,
              binary_mod_pb=seasonal_mods_pb$binary_mod,
              qt_mod_pb=seasonal_mods_pb$qt_mod))
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

KM_gg <- function(survfit_obj,title=''){
  survfit_dat <- tibble(time=rep(survfit_obj$time,survfit_obj$n.event),
                        surv=rep(survfit_obj$surv,survfit_obj$n.event),
                        strata=rep(paste0(c(0,1,2),' (n=',survfit_obj$n,')'),survfit_obj$n)) %>%
    distinct()
  ggplot(survfit_dat,aes(time,surv,col=factor(strata))) +
    geom_step() +
    scale_y_continuous(limits=c(0,1),expand=c(0,0)) +
    scale_color_manual(values=c("#BC3C29FF","#0072B5FF","#E18727FF"),name='Genotype') +
    ggtitle(title) +
    xlab('Time in months') +
    ylab('Survival') +
    theme_classic() +
    theme(legend.position=c(0.8,0.8))
}
