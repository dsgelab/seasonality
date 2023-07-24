library(readr)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
source('utils.R')

get_monte_carlo_PTR_CI <- function(mod,nr_iter,adjustment,seasonal_spline_avg){
  relevant <- grep('EVENT_MONTH',names(mod$coefficients))
  V_f <- vcov(mod)[relevant,relevant]
  s_hat <- coef(mod)[relevant]
  prediction_grid <- data.frame(EVENT_MONTH_DEC=seq(1,12.99,by=0.01))
  lp <- predict(mod, newdata = prediction_grid , type = "lpmatrix")[,relevant]
  s_sample <- mvrnorm(n=nr_iter,mu=s_hat,Sigma=V_f)
  smooth_sample <- lp %*% t(s_sample)
  metrics <- apply(smooth_sample,2,function(x){
    c('ptr'=exp(max(x))/exp(min(x)),
      'peak_month'= prediction_grid$EVENT_MONTH_DEC[which.max(x)],
      'trough_month'= prediction_grid$EVENT_MONTH_DEC[which.min(x)])
  })
  ptr_CI <- quantile(metrics['ptr',],probs=c(0.025,0.975))
  month_CI <- apply(metrics[c('peak_month','trough_month'),],1,function(x){
    #detect if vector breaks year cycle (lateest third and first third of the year included)
    if(any(x<4) & any(x>8)){
      x <- ifelse(x>6,x-12,x)
    }
    q_x <- quantile(x,probs=c(0.025,0.975))
    q_x <- ifelse(q_x<1,q_x+12,q_x)
    return(q_x)
  })
  return(list(ptr=ptr_CI,month=month_CI))
}

summarise_seasonality <- function(mod){
  seasonal_grid <- data.frame(EVENT_MONTH_DEC=seq(1,12.99,by=0.01))
  seasonal_component_grid <- predict(mod,newdata=seasonal_grid,type='terms')
  seasonal_component_grid_est <- seasonal_component_grid[,'s(EVENT_MONTH_DEC)']
  idx_peak <- which.max(seasonal_component_grid_est)
  val_peak <- exp(seasonal_component_grid_est[idx_peak])
  idx_trough <- which.min(seasonal_component_grid_est)
  val_trough <- exp(seasonal_component_grid_est[idx_trough])
  month_peak <- seasonal_grid$EVENT_MONTH_DEC[idx_peak]
  month_trough <- seasonal_grid$EVENT_MONTH_DEC[idx_trough]
  ptr_est <- val_peak / val_trough
  CI <- get_monte_carlo_PTR_CI(mod=mod,nr_iter=100)
  tibble(ptr_est=ptr_est,
         ptr_lower=CI$ptr['2.5%'],
         ptr_upper=CI$ptr['97.5%'],
         smooth_pval=summary(mod)$s.table['s(EVENT_MONTH_DEC)','p-value'],
         month_peak=month_peak,
         month_peak_lower=CI$month['2.5%','peak_month'],
         month_peak_upper=CI$month['97.5%','peak_month'],
         month_trough=month_trough,
         month_trough_lower=CI$month['2.5%','trough_month'],
         month_trough_upper=CI$month['97.5%','trough_month'],
         val_peak=val_peak,
         val_trough=val_trough)
}

diagnosis_dates <- read_tsv('data/FINNGEN_endpoints_individual_dates.txt') %>%
                   filter(EVENT_YEAR>=1998,EVENT_YEAR<=2019)
seasonal_smooth_terms <- read_tsv('data/FINREGISTRY_seasonal_splines.txt')

intervene_path <- '/finngen/red/Zhiyu/Share/R10_PRS/'
endpoint_ids <- c('Asthma'='ASTHMA_NAS','I9_STR-wSampleOverlap'='I9_STR','G6_AD_WIDE'='G6_AD_WIDE','I9_AF'='I9_AF','N14_CHRONKIDNEYDIS'='N14_CHRONKIDNEYDIS','Rheumatoid_Arthritis'='M13_RHEUMA','T1D'='T1D',
                  'I9_CHD'='I9_CHD','Inflammatory_Bowel_Disease'='K11_IBD_STRICT','Osteoporosis'='M13_OSTEOPOROSIS','Sleep_Apnoea'='G6_SLEEPAPNO','T2D'='T2D',
                  'Appendicitis'='K11_APPENDACUT','GOUT'='GOUT_STRICT','I9_HEARTFAIL_NS'='I9_HEARTFAIL_ALLCAUSE','MDD'='F5_DEPRESSIO','Stroke'='C_STROKE')
# endpoint <- 'F5_DEPRESSIO'
# endpoint_PRS_path <- '/finngen/library-red/finngen_R11/prs_1.0/data/finngen_R11_PGC_UKB_depression_genome-wide.txt.sscore'

PRS_seasonality_results <- lapply(1:length(endpoint_ids),function(i){
                            diagnosis_dates_endpoint <- filter(diagnosis_dates,ENDPOINT==endpoint_ids[i])
                            PRS_dat <- read_tsv(paste0(intervene_path,names(endpoint_ids)[i],'.sscore')) %>%
                                       mutate(SCORE1=scale(SCORE1_AVG))

                            seasonal_smooth_term_endpoint <- filter(seasonal_smooth_terms,ENDPOINT==endpoint_ids[i]) %>%
                                                    mutate(month=ifelse(month<1,month+12,month)) %>%
                                                    distinct(month,.keep_all=T)

                            PRS_time_dat <- inner_join(diagnosis_dates_endpoint,PRS_dat,by=c('FINNGENID'='#FID')) %>%
                              mutate(EVENT_MONTH_DEC=get_month_dec(EVENT_DATE,year_start=1998,year_end=2019),
                                     seasonal_val=approx(x=seasonal_smooth_term_endpoint$month,
                                                         y=seasonal_smooth_term_endpoint$seasonal_val,
                                                         xout=EVENT_MONTH_DEC)$y,
                                     season=ifelse(seasonal_val>0,'high','low'))

                            PRS_gam <- gam(SCORE1~s(EVENT_MONTH_DEC,bs = 'cc',k=6),data=PRS_time_dat,family=gaussian())
                            PRS_season <- summary(lm(SCORE1~season,data=PRS_time_dat))

                            summarise_seasonality(PRS_gam) %>%
                            mutate(ENDPOINT=endpoint_ids[i]) %>%
                            mutate(season_effect=PRS_season$coefficients['seasonlow',1],
                                   season_pval=PRS_season$coefficients['seasonlow',4]) %>%
                            dplyr::select(ENDPOINT,season_pval,everything())
                          }) %>% bind_rows()

mutate(PRS_time,EVENT_MONTH=floor(EVENT_MONTH_DEC)) %>%
  ggplot(aes(x=as.factor(EVENT_MONTH),y=SCORE1)) +
  geom_boxplot()

plot(PRS_gam)

