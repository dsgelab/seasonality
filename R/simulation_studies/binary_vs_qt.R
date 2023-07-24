library(readr)
library(data.table)
library(gratia)
source('seasonality_dashboard/utils.R')

GWAS_info_dat <- read_tsv('data/finngen_R10_finngen_R10_analysis_data_finngen_R10_pheno_n.tsv') %>%
                 filter(num_gw_significant>1)

monthly_counts_FinRegistry  <- read_tsv('data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
                               filter_monthly_counts(.,GWAS_info_dat)
seasonal_spline <- read_tsv('data/FINREGISTRY_seasonal_splines.txt')

endpoint_id='F5_DEPRESSIO'
dat_endpoint <- filter(monthly_counts_FinRegistry,ENDPOINT==endpoint_id)
seasonal_spline_endpoint <- filter(seasonal_spline,ENDPOINT==endpoint_id)
mod_list <- run_seasonality_gam(dat_endpoint)


n_cases <- filter(GWAS_info_dat,phenocode==endpoint_id) %>% .$num_cases
MAF <- 0.5
beta <- 0.5
x <- rbinom(n_cases,size=2,prob=MAF)
y <- rnorm(n_cases,mean=beta*x,sd=1)



count_dat <- as_tibble(simulate(mod_list$seasonal)) %>%
             rename(count=V1) %>%
             bind_cols(expand_grid(year=seq(1998,2019),
                                   month=seq_len(12)),.) %>%
             mutate(date=ym(paste0(year,'-',month)))

all_dates <- seq(ymd('1998-01-01'),ymd('2019-12-31'),by='1 day')
month_days <- tibble(year=year(all_dates),month=month(all_dates),day=day(all_dates)) %>%
  group_by(year,month) %>%
  summarise(num_days=max(day))

date_dat <- lapply(1:nrow(count_dat),function(i){
              if(i<nrow(count_dat)){
                dates <- head(seq(count_dat$date[i],count_dat$date[i+1],by='1 day'),-1) #remove first day from next month
              }else{
                dates <- head(seq(count_dat$date[i],ym('2020-01'),by='1 day'),-1)
              }
              sim_dates <- sample(dates,size=count_dat$count[i],replace=T)
              tibble(year=count_dat$year[i],
                     month=count_dat$month[i],
                     day=day(sim_dates),
                     date=sim_dates)
            }) %>% bind_rows() %>%
            sample_n(n_cases) %>%
            inner_join(month_days,by=c('year','month')) %>%
            mutate(event_month_dec=month + day/num_days) %>%
            mutate(event_month_dec=ifelse(event_month_dec>12.5,
                                          event_month_dec-12,
                                          event_month_dec))
eps <- 1e-6
pheno_dat_endpoint <- mutate(date_dat,
                             seasonal_val=approx(x=seasonal_spline_endpoint$month,
                                                 y=seasonal_spline_endpoint$seasonal_val,
                                                 xout=event_month_dec)$y,
                             seasonal_val_beta=(seasonal_val-min(seasonal_val))/(max(seasonal_val)-min(seasonal_val)),
                             seasonal_val_beta=case_when(abs(seasonal_val_beta-0)<eps ~ seasonal_val_beta+eps,
                                                        abs(seasonal_val_beta-1)<eps ~ seasonal_val_beta-eps,
                                                        TRUE ~ seasonal_val_beta),
                             seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                             seasonal_val_binary=as.integer(seasonal_val > 0))

phi=1.068
mu=0.5848
a=mu*phi
b=phi*(1-mu)
plot(dbeta(seq(0,1,by=0.001),shape1=a,shape2=b),type='l')

x <- 1.4
gamma <- x*mu/(1-x*mu)
a1 <- (a+b)*gamma/(1+gamma)
b1 <- a+b-a1

plot(dbeta(seq(0,1,by=0.001),shape1=a,shape2=b),type='l')
lines(dbeta(seq(0,1,by=0.001),shape1=a1,shape2=b1),type='l')


ggplot(filter(dat,year<2005),aes(date,count)) +
geom_line()


library(flexsurv)
hseasonal <- function(x,g,beta,phi,s,P){
  beta[1] + beta[2]*g + beta[3]*sin(2*pi*(x-s) / P + phi) + beta[4]*g*sin(2*pi*(x-s)/P + phi)
}

log_hseasonal <- function(x,g,beta,phi,s,P){
  beta[1] + beta[2]*g + beta[3]*sin(2*pi*(x-s) / P + phi) + beta[4]*g*sin(2*pi*(x-s)/P + phi)
}

#Integrand of hazard function
I_seasonal <- function(x,g,beta,phi,s,P){
  (beta[1] + beta[2]*g)*(x-s) - (P*(beta[3]+beta[4]*g)/(2*pi))*cos(2*pi*(x-s) / P + phi)
}

Hseasonal <- function(x,g,beta,phi,s,P){
  I_seasonal(x,g,beta,phi,s,P) - I_seasonal(s,g,beta,phi,s,P)
}

P <- 12
s <- 0.5
t <- seq(s,P+s,by=0.01)
#beta=c(0.1,0,0,-0.025)
beta=c(-1.5,0,0,0.25)
phi <- pi/2
gt <- 0:2
dat_list <- lapply(1:length(gt),function(i){
  haz <- exp(hseasonal(t,gt[i],beta,phi,s,P))
  cum_haz <- Hseasonal(t,gt[i],beta,phi,s,P)
  cum_haz2 <- cumsum(haz*0.01)
  surv <- exp(-cum_haz2)
  tibble(t=t,gt=gt[i],haz=haz,cum_haz=cum_haz,cum_haz2=cum_haz2,surv=surv)
})

ggplot(bind_rows(dat_list),aes(t,surv,col=as.factor(gt))) +
geom_line() +
scale_y_continuous(limits=c(0,1)) +
theme_classic()

ggplot(dat,aes(t,cum_haz2,col=gt)) +
geom_line() +
theme_classic()

ggplot(bind_rows(dat_list),aes(t,haz,col=as.factor(gt))) +
  geom_line() +
  theme_classic()

n <- 40000
MAF <- 0.05
gt <- rbinom(n,size=2,prob=MAF)
u <- runif(n = n,min = 0,max=1-min(dat$surv))
event_month_dec <- vector('numeric',n)
for(i in 1:n){
  event_month_dec[i] <- dat_list[[gt[i]+1]]$t[which.min(abs(dat_list[[gt[i]+1]]$surv-u[i]))]
}
sim_dat <- tibble(gt=gt,event_month_dec=event_month_dec,surv=u)


KM_gg <- function(survfit_obj,title=''){
  survfit_dat <- tibble(time=rep(survfit_obj$time,survfit_obj$n.event),
                        surv=rep(survfit_obj$surv,survfit_obj$n.event),
                        strata=rep(c(0,1,2),survfit_obj$n)) %>%
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

surv_obj_sim <- Surv(time=sim_dat$event_month_dec,event=rep(1,n))
survfit_sim <- survfit(surv_obj_sim ~ gt)
KM_gg(survfit_sim)

lambda_vec <- seq(0,3*28)
diagnosis_lag <- rpois(n,lambda=tail(lambda_vec,1))/28
perturb <- sample(seq_len(15),size = n,replace=T)/28
perturb <- ifelse(runif(n)<0.5,-1*perturb,perturb)
sim_dat$event_month_dec_pb <- sim_dat$event_month_dec+diagnosis_lag+perturb
sim_dat$event_month_dec_pb <- ifelse(sim_dat$event_month_dec_pb > 12.5,sim_dat$event_month_dec_pb-12,
                              ifelse(sim_dat$event_month_dec_pb < 0.5,sim_dat$event_month_dec_pb+12,sim_dat$event_month_dec_pb))

surv_obj_sim_pb <- Surv(time=sim_dat$event_month_dec_pb,event=rep(1,n))
survfit_sim_pb <- survfit(surv_obj_sim_pb ~ gt)
KM_gg(survfit_sim_pb)

KM_gg(survfit_sim,'True disease onset') + (KM_gg(survfit_sim_pb,title='True disease onset with perturbation') + theme(legend.position='none'))


monthly_counts <- mutate(sim_dat,
                                 month=ifelse(event_month_dec<1,12,floor(event_month_dec)),
                                 month_pb=ifelse(event_month_dec_pb<1,12,floor(event_month_dec_pb))) %>%
                  pivot_longer(cols=c('month','month_pb')) %>%
                  group_by(name,value) %>%
                  summarise(count=length(value))
