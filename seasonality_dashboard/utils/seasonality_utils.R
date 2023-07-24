library(MASS,exclude = "select")
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(mgcv)
library(readxl)

#### Data processing ####
get_endpoint_mapping <- function(path,seasonal_summary){
  endpoint_groups <- read_excel(path,sheet=2) %>%
                      mutate(code=gsub('#','',code)) %>%
                      pivot_longer(cols=c('CHAPTER','OTHER'),names_to='descr_type',values_to='descr',values_drop_na = T)
  endpoint_mapping <- mutate(seasonal_summary,prefix=gsub('_.*$','',ENDPOINT)) %>%
                      select(ENDPOINT,prefix) %>%
                      inner_join(endpoint_groups,by=c('prefix'='code'))
  return(endpoint_mapping)
}

filter_monthly_counts <- function(monthly_counts,endpoint_mapping_path){
  monthly_counts_filtered <- filter(monthly_counts,EVENT_YEAR>=1998 & EVENT_YEAR<=2019)
  month_year_endpoint_template <- expand_grid(ENDPOINT=unique(monthly_counts_filtered$ENDPOINT),
                                              EVENT_YEAR=seq(min(monthly_counts_filtered$EVENT_YEAR),
                                                             max(monthly_counts_filtered$EVENT_YEAR)),
                                              EVENT_MONTH=seq_len(12))
  #endpoint_mapping <- get_endpoint_mapping(path=endpoint_mapping_path,seasonal_summary=distinct(monthly_counts_filtered,ENDPOINT))
  monthly_counts_filtered <- left_join(month_year_endpoint_template,
                                       monthly_counts_filtered,
                                       by=c('ENDPOINT','EVENT_YEAR','EVENT_MONTH')) %>%
                              mutate(COUNT=ifelse(is.na(COUNT),0L,COUNT),
                                     EVENT_DATE=ym(paste0(EVENT_YEAR,'-',EVENT_MONTH))) %>%
                              #left_join(select(endpoint_mapping,ENDPOINT,prefix,brief),by='ENDPOINT') %>%
                              #filter(is.na(prefix) | !(prefix %in% c('OTHER','VWXY20'))) %>%
                              arrange(ENDPOINT,EVENT_YEAR,EVENT_MONTH)
  #all endpoints have to have more than 0 diagnoses each year.
  monthly_counts_filtered  <- group_by(monthly_counts_filtered,ENDPOINT,EVENT_YEAR) %>%
                              summarise(count=sum(COUNT)) %>%
                              group_by(ENDPOINT) %>%
                              filter(all(count>0)) %>%
                              distinct(ENDPOINT) %>%
                              inner_join(monthly_counts_filtered,.,by='ENDPOINT')
  return(monthly_counts_filtered)
}

#### Fit seasonality GAM #####
get_seasonal_spline_adj <- function(seasonal_spline_avg,months){
  approx(x=seasonal_spline_avg$month,
         y=seasonal_spline_avg$avg_seasonal_val,
         xout=months)$y
}

get_grid <- function(type,adjustment='',by_year=0.1,by_month=0.1,seasonal_spline_avg=NULL){
  if(type=='seasonal'){
   grid_dat <-data.frame(EVENT_YEAR=1998, EVENT_MONTH=seq(0.5,12.5,by=by_month))
  }else if(type=='trend'){
   grid_dat <- data.frame(EVENT_YEAR=seq(1998,2019,by=by_year),EVENT_MONTH=1)
  }else{
    stop('Type not recognized')
  }
  grid_dat$nr_days <- 30
  if(grepl('mean',adjustment)){
    grid_dat <- mutate(grid_dat,avg_seasonal_val=get_seasonal_spline_adj(seasonal_spline_avg,EVENT_MONTH))
  }else if(adjustment=='binary'){
    grid_dat <- mutate(grid_dat,july=as.integer(floor(EVENT_MONTH)==7),december=as.integer(floor(EVENT_MONTH)==12 | EVENT_MONTH<1))

  }
  return(grid_dat)
}

run_seasonality_gam <- function(dat,a,b,adjustment='',adj_month_length=F,mod_null='adj',seasonal_spline_type='cp',k_seasonal=6,seasonal_spline_avg=NULL){
  dates <- seq(ymd('1998-01-01'),ymd('2019-12-31'),by='1 day')
  dates_dat <- tibble(EVENT_YEAR=year(dates),EVENT_MONTH=month(dates),EVENT_DAY=day(dates)) %>%
               group_by(EVENT_YEAR,EVENT_MONTH) %>%
               summarise(nr_days=max(EVENT_DAY),.groups = 'drop')
  k_trend <- 6
  if(adj_month_length){
    dat <- inner_join(dat,dates_dat,by=c('EVENT_YEAR','EVENT_MONTH'))
  }else{
    dat$nr_days <- 30
  }
  offset_null <- 'log(nr_days)'
  if(adjustment=='mean_covariate'){
    dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,dat$EVENT_MONTH)
    offset_adj <- ''
    cov_adj <- ' + avg_seasonal_val'
  }else if(adjustment=='binary'){
    dat$july <- as.integer(floor(dat$EVENT_MONTH)==7)
    dat$december <- as.integer(floor(dat$EVENT_MONTH)==12)
    offset_adj <- ''
    cov_adj <- ' + july + december'
  }else if(adjustment=='mean_offset'){
    dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,dat$EVENT_MONTH)
    offset_adj <- ' + avg_seasonal_val'
    cov_adj <- ''
  }else{
    offset_adj <- ''
    cov_adj <- ''
  }
  f_start <- 'COUNT ~ s(EVENT_YEAR, k=k_trend, bs="ps")'
  if(mod_null=='adj'){
    f_null <- paste0(f_start,' + offset(',offset_null,offset_adj,')',cov_adj)
    f_seasonal <- f_null
  }else if(mod_null=='noadj'){
    f_null <- paste0(f_start,' + offset(',offset_null,')')
    f_seasonal <- paste0(f_start,' + offset(',offset_null,offset_adj,')',cov_adj)
  }else{
    stop('mod_null argument not recognized')
  }
  f_seasonal <- paste0(f_seasonal,' + s(EVENT_MONTH, k=k_seasonal, bs="',seasonal_spline_type,'")')
  family <- quasipoisson()
  #family <- nb()
  seasonal <- gam(as.formula(f_seasonal), family=family, data=dat, knots=list(EVENT_MONTH = c(a-0.5, b-0.5)), scale=-1,method = "REML")
  null <- gam(as.formula(f_null),family=family, data=dat, scale=-1,method = "REML")
  return(list(seasonal=seasonal, null=null))
}

summarise_seasonality <- function(dat,a,b,adjustment='',adj_month_length=F,mod_null='adj',seasonal_spline_type='cp',k_seasonal=6,seasonal_spline_avg=NULL){
  mod_list <- run_seasonality_gam(dat=dat,a=a,b=b,adjustment=adjustment,adj_month_length=adj_month_length,mod_null=mod_null,seasonal_spline_type=seasonal_spline_type,
                                  k_seasonal=k_seasonal,seasonal_spline_avg=seasonal_spline_avg)
  seasonal_grid <- get_grid(type='seasonal',adjustment=adjustment,seasonal_spline_avg=seasonal_spline_avg,by_month = 0.01)
  seasonal_component_monthly <- predict(mod_list$seasonal,
                                        newdata=seasonal_grid,
                                        type='terms')
  seasonal_component_monthly_est <- seasonal_component_monthly[,'s(EVENT_MONTH)']
  idx_peak <- which.max(seasonal_component_monthly_est)
  val_peak <- exp(seasonal_component_monthly_est[idx_peak])
  idx_trough <- which.min(seasonal_component_monthly_est)
  val_trough <- exp(seasonal_component_monthly_est[idx_trough])
  month_peak <- seasonal_grid$EVENT_MONTH[idx_peak]
  month_trough <- seasonal_grid$EVENT_MONTH[idx_trough]
  ptr_est <- val_peak / val_trough
  CI <- get_monte_carlo_PTR_CI(mod=mod_list$seasonal,nr_iter=100,adjustment=adjustment,seasonal_spline_avg = seasonal_spline_avg)
  tibble(ptr_est=ptr_est,
         ptr_lower=CI$ptr['2.5%'],
         ptr_upper=CI$ptr['97.5%'],
         smooth_pval=summary(mod_list$seasonal)$s.table['s(EVENT_MONTH)','p-value'],
         month_peak=month_peak,
         month_peak_lower=CI$month['2.5%','peak_month'],
         month_peak_upper=CI$month['97.5%','peak_month'],
         month_trough=month_trough,
         month_trough_lower=CI$month['2.5%','trough_month'],
         month_trough_upper=CI$month['97.5%','trough_month'],
         val_peak=val_peak,
         val_trough=val_trough,
         deviance=mod_list$seasonal$deviance,
         deviance_expl=(mod_list$null$deviance-mod_list$seasonal$deviance)/mod_list$null$deviance,
         dispersion=mod_list$seasonal$scale)
}

extract_seasonal_spline <- function(dat,a,b,adjustment='',adj_month_length=F,seasonal_spline_type='cp',k_seasonal=6,seasonal_spline_avg=NULL){
  mod_list <- run_seasonality_gam(dat,a=a,b=b,adjustment=adjustment,adj_month_length=adj_month_length,seasonal_spline_type=seasonal_spline_type,k_seasonal=k_seasonal,seasonal_spline_avg = seasonal_spline_avg)
  months_vec <- seq(a-0.5,b-0.5,by=0.01)
  seasonal_pred <- predict(mod_list$seasonal,newdata=get_grid(type='seasonal',adjustment=adjustment,by_month = 0.01,seasonal_spline_avg=seasonal_spline_avg),type='terms')[,'s(EVENT_MONTH)']
  return(tibble(month=months_vec,seasonal_val=seasonal_pred))
}

get_monte_carlo_PTR_CI <- function(mod,nr_iter,adjustment,seasonal_spline_avg){
  relevant <- grep('EVENT_MONTH',names(mod$coefficients))
  V_f <- vcov(mod)[relevant,relevant]
  s_hat <- coef(mod)[relevant]
  prediction_grid <- get_grid(type='seasonal',by_month = 0.01,adjustment=adjustment,seasonal_spline_avg=seasonal_spline_avg)
  lp <- predict(mod, newdata = prediction_grid , type = "lpmatrix")[,relevant]
  s_sample <- mvrnorm(n=nr_iter,mu=s_hat,Sigma=V_f)
  smooth_sample <- lp %*% t(s_sample)
  month_grid_adjusted <- ifelse(prediction_grid$EVENT_MONTH<1,prediction_grid$EVENT_MONTH+12,prediction_grid$EVENT_MONTH)
  metrics <- apply(smooth_sample,2,function(x){
    c('ptr'=exp(max(x))/exp(min(x)),
      'peak_month'= month_grid_adjusted[which.max(x)],
      'trough_month'= month_grid_adjusted[which.min(x)])
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

#### Seasonality plots #####
theme_custom <- theme(plot.title = element_text(size=16),
                      axis.title.x = element_text(size=16),
                      axis.text.x = element_text(size=14),
                      axis.title.y = element_text(size=16),
                      axis.text.y = element_text(size=14),
                      legend.title = element_text(size=16),
                      legend.text = element_text(size=14))

compare_gam_fits <- function(mod_list,monthly_counts,endpoint,adjustment=''){
  years <- unique(monthly_counts$EVENT_YEAR)
  season_dat <- bind_rows(list(summer=tibble(start=ym(paste0(years,'-',4)),
                                             end=ym(paste0(years,'-',8))),
                               winter=tibble(start=ym(paste0(years,'-',8)),
                                             end=ym(paste0(years+1,'-',4)))),.id='season')
  season_dat <- bind_rows(season_dat,tibble(season='winter',
                                            start=ym(paste0(years[1],'-',1)),
                                            end=ym(paste0(years[1],'-',4))))
  monthly_counts$null_pred <-predict(mod_list$null)
  monthly_counts$seasonal_pred <-predict(mod_list$seasonal)
  monthly_counts <- mutate(monthly_counts,month=1:n())
  p <- ggplot(monthly_counts) +
        geom_point(aes(x=EVENT_DATE,
                       y=COUNT)) +
        geom_line(aes(x=EVENT_DATE,
                      y=exp(seasonal_pred))) +
          geom_line(aes(x=EVENT_DATE,
                        y=exp(null_pred)),
                    col='red') +
        geom_rect(data=season_dat,aes(xmin=start,
                                      xmax=end,
                                      ymin=min(monthly_counts$COUNT),
                                      ymax=max(monthly_counts$COUNT),
                                      fill=season),alpha=0.1) +
        scale_y_continuous(expand=c(0.01,0)) +
        scale_x_date(expand=c(0.005,0),date_breaks='2 years',date_labels='%Y') +
        scale_fill_manual(values=c('gold3','turquoise')) +
        xlab(paste0('First diagnosis of: ',endpoint)) +
        ylab('Count') +
        theme_classic() +
        theme(legend.position='none',
              axis.title.x=element_text(size=14),
              axis.text.x=element_text(size=14),
              axis.title.y=element_text(size=14),
              axis.text.y=element_text(size=14))
  return(p)
}

trend_plot <- function(dat,mod_list,adjustment,seasonal_spline_avg=NULL){
  years_vec <- seq(1998,2019,by=0.1)
  trend_pred <- predict(mod_list$seasonal,newdata=get_grid(type='trend',adjustment=adjustment,seasonal_spline_avg=seasonal_spline_avg),type='terms',unconditional=T,se.fit=T)
  trend_pred_dat <- tibble(year=years_vec,
                           est=trend_pred$fit[,'s(EVENT_YEAR)'],
                           lower=est - 1.96*trend_pred$se.fit[,'s(EVENT_YEAR)'],
                           upper=est + 1.96*trend_pred$se.fit[,'s(EVENT_YEAR)'])
  ggplot(trend_pred_dat) +
  geom_line(aes(x=year,
                y=est)) +
  geom_ribbon(aes(x=year,ymin=lower,ymax=upper),alpha=0.3,fill='blue') +
  geom_hline(yintercept=0,linetype='dashed',color='red') +
  scale_x_continuous(limits = c(1998,2020),expand=c(0.03,0)) +
  xlab('Year') +
  ylab('Smooth term') +
  ggtitle('Trend component') +
  theme_bw() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14))
}

seasonality_plot <- function(dat,mod_list,a,b,adjustment,seasonal_spline_avg=NULL){
  months_vec <- seq(a-0.5,b-0.5,by=0.1)
  seasonal_pred <- predict(mod_list$seasonal,newdata=get_grid(type='seasonal',adjustment=adjustment,seasonal_spline_avg=seasonal_spline_avg),type='terms',unconditional=T,se.fit=T)
  seasonal_pred_dat <- tibble(month=months_vec,
                              est=seasonal_pred$fit[,'s(EVENT_MONTH)'],
                              lower=est - 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'],
                              upper=est + 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'])
  ggplot(seasonal_pred_dat) +
  geom_line(aes(x=month,
                y=est)) +
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.3,fill='blue') +
  geom_hline(yintercept=0,linetype='dashed',color='red') +
  scale_x_continuous(breaks=seq(a,b-1)) +
  xlab('Month number') +
  ylab('Smooth term') +
  ggtitle('Seasonal component') +
  theme_bw() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14))
}
