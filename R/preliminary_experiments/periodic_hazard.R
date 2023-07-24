library(survival)
library(dplyr)
library(ggplot2)
library(flexsurv)
N <- 1000
year <- sample(0:5,size = N,replace=T)
day_trend <- 1.1+cos(2*pi*(1:365)/365)
day_distr <- day_trend/sum(day_trend)
day <- sample(1:365,size=N,prob=day_distr,replace = T)
times <- year*365 + day
plot(survfit(Surv(times,event=rep(1,N))~1))

#
#log(lambda) <- alpha+beta*cos(2*pi*t/365)

log_hazard <- function(t,g,beta,P=365){
  beta[1] + beta[2]*g +
  beta[3]*cos(2*pi*t/P) +  beta[4]*sin(2*pi*t/P) +
  beta[5]*g*cos(2*pi*t/P) +  beta[6]*g*sin(2*pi*t/P)
}

N_days <- 365
t <- seq(0,10*N_days)
beta=c(-7.5,0,-1,1,0.5,0)
# alpha <- -3.8
# beta <- 2
# gamma <- 0
# delta <- 0.2
gt <- 0:2
phi <- 0
dat <- lapply(1:length(gt),function(i){
        log_haz <- log_hazard(t,gt[i],beta)
        #log_haz <- exp(alpha+gamma*gt[i])*exp((beta+delta*gt[i])*cos(phi+2*pi*t/N_days))
        cum_haz <- cumsum(exp(log_haz))
        surv <- exp(-cum_haz)
        tibble(t=t,gt=gt[i],log_haz=log_haz,haz=exp(log_haz),cum_haz=cum_haz,surv=surv)
      }) %>% bind_rows() %>% mutate(gt=factor(gt))

ggplot(dat,aes(t,surv,col=gt)) +
geom_line() +
theme_classic()

ggplot(dat,aes(t,cum_haz,col=gt)) +
  geom_line() +
  theme_classic()

ggplot(dat,aes(t,haz,col=gt)) +
  geom_line() +
  theme_classic()

ggplot(dat,aes(t,log(haz),col=gt)) +
  geom_line() +
  theme_classic()


MAF <- 0.3
n <- 2000
gt <- sample(0:2,size=n,prob = c(1-MAF,MAF*(1-MAF),MAF^2),replace=T)
times <- sapply(1:n,function(i){
  dat_iter <- dat[dat$gt==gt[i],]
  dat_iter$t[which.min(abs(dat_iter$surv-runif(1)))]
})
plot(survfit(Surv(times,event=rep(1,n)) ~ gt))


c1 <- -1
c2 <- 1
t <- seq(0,10,by=0.01)
plot(t,c1*cos(t) + c2*sin(t))


