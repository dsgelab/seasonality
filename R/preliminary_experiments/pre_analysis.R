library(dplyr)
library(mgcv)
library(ggplot2)
library(brms)
source('utils.R')
# Assumptions
N_days <- 365
t <- 1:N_days
dates <- seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by="+1 day")
mod_type <- ''
dat <- tibble(t=t,date=dates)
# Seasonal variation of disease diagnoses
A <- 10
shift <- 100
baseline <- 30
y_mu_count <- baseline + A*cos(2*pi*t/N_days)
y_count <- rnbinom(N_days,mu=y_mu_count,size=50)
dat_count <- bind_cols(dat,tibble(y=y_count,y_mu=y_mu_count))
if(mod_type=='brms'){
  mod_count <- brm(y ~ s(t,bs='cc'),family=negbinomial(),data=dat)
}else{
  mod_count <- gam(y ~ s(t, bs = "cc"),family=nb(),data=dat_count)
}

dat_count <- bind_cols(dat_count,get_pred_gamm(mod_count,dat_count,type=mod_type,tr=exp))


ggplot(data=dat_count) +
geom_point(aes(date,y),size=3,alpha=0.5) +
geom_line(aes(date,mean,col='Predicted'),size=2) +
geom_ribbon(aes(date,ymin=lower,ymax=upper,col='Predicted'),alpha=0.1) +
geom_line(aes(date,y_mu,col='Actual'),size=2) +
scale_color_manual(name='',
                   breaks=c('Actual', 'Predicted'),
                   values=c('Actual'='#BC3C29FF', 'Predicted'='#0072B5FF')) +
xlab('Date') +
ylab('Diagnosis count') +
theme_classic() +
theme_font



### Seasonal distribution of PRS
y_mu_prs <- 2*sin((2*pi*(t-shift))/N_days)
y_prs <- rnorm(N_days,mean=y_mu_prs,sd=1)
dat_prs <- bind_cols(dat,tibble(y=y_prs,y_mu=y_mu_prs))

#We are not able to extract seasonal prs effect through gam modeling alone
if(mod_type=='brms'){
  mod_prs <- brm(y ~ s(t,bs='cc'),family=gaussian(),data=dat_prs)
}else{
  mod_prs <- gam(y ~ s(t,bs='cc'),family=gaussian(),data=dat_prs)
}
dat_prs <- bind_cols(dat_prs,get_pred_gamm(mod_prs,dat_prs,type=mod_type))

ggplot(data=dat_prs) +
geom_point(aes(date,y),size=3,alpha=0.5) +
geom_line(aes(date,mean,col='Predicted'),size=2) +
geom_ribbon(aes(date,ymin=lower,ymax=upper,col='Predicted'),alpha=0.1) +
geom_line(aes(date,y_mu,col='Actual'),size=2) +
scale_color_manual(name='',
                   breaks=c('Actual', 'Predicted'),
                   values=c('Actual'='#BC3C29FF', 'Predicted'='#0072B5FF')) +
xlab('Date') +
ylab('Polygenic risk score') +
theme_classic() +
theme_font


# Individual level data
gt <- sample(0:2,prob = c(0.8,0.16,0.04),10000,replace=T)
diag <- vector('numeric',length(gt))
diag[gt==0] <- ceiling(365*rbeta(sum(gt==0),shape1 = 0.6,shape2=0.6))
diag[gt==1] <- ceiling(365*rbeta(sum(gt==1),shape1 = 0.8,shape2=0.8))
diag[gt==2] <- ceiling(365*rbeta(sum(gt==2),shape1 = 1,shape2=1))

plot(dat$t,table(diag))


# if(gt==0){
#
# }
# prs <- rnorm(1000,mean=0,sd=1)
# prs_quant <- pnorm(prs)
#
# sapply(1:length(prs),function(i){
#   pers_trend <- cos(2*pi*t/N_days) + baseline
#   pers_trend <- pers_trend/sum(pers_trend)
#
# })
# trend_normalized <- trend/sum(trend)
library(lubridate)
counts <- c(2500,2450,2400,23)
seq(ymd('1997-01-01'),ymd('2001-12-01'),by="month")
dat <- tibble(date=)
