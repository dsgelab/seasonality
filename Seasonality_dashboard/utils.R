filter_monthly_counts <- function(monthly_counts,GWAS_info_dat){
  monthly_counts_filtered <- filter(monthly_counts,EVENT_YEAR>=1998 & EVENT_YEAR<=2019)
  month_year_endpoint_template <- expand_grid(ENDPOINT=unique(monthly_counts_filtered$ENDPOINT),
                                              EVENT_YEAR=seq(min(monthly_counts_filtered$EVENT_YEAR),
                                                             max(monthly_counts_filtered$EVENT_YEAR)),
                                              EVENT_MONTH=seq_len(12))
  monthly_counts_filtered <- left_join(month_year_endpoint_template,
                                       monthly_counts_filtered,
                                       by=c('ENDPOINT','EVENT_YEAR','EVENT_MONTH')) %>%
    mutate(COUNT=ifelse(is.na(COUNT),0L,COUNT),
           EVENT_DATE=ym(paste0(EVENT_YEAR,'-',EVENT_MONTH))) %>%
    inner_join(select(GWAS_info_dat,ENDPOINT=phenocode,num_gw_significant),by='ENDPOINT') %>%
    arrange(ENDPOINT,EVENT_YEAR,EVENT_MONTH)
  return(monthly_counts_filtered)
}

get_seasonal_spline_adj <- function(seasonal_spline_avg,months){
  approx(x=seasonal_spline_avg$month,
         y=seasonal_spline_avg$avg_seasonal_val,
         xout=months)$y
}

get_grid <- function(dat,type,seasonal_spline_avg=NULL){
  if(type=='seasonal'){
   grid_dat <-data.frame(EVENT_YEAR=1998, EVENT_MONTH=seq(0.5,12.5,by=0.1))
  }else if(type=='trend'){
   grid_dat <- data.frame(EVENT_YEAR=seq(1998,2019,by=0.1),EVENT_MONTH=1)
  }else{
    stop('Type not recognized')
  }
  grid_dat$nr_days <- 30
  if(!is.null(seasonal_spline_avg)){
    grid_dat <- mutate(grid_dat,avg_seasonal_val=get_seasonal_spline_adj(seasonal_spline_avg,EVENT_MONTH))
  }
  return(grid_dat)
}

run_seasonality_gam <- function(dat,a,b,seasonal_spline_avg=NULL){
  k_trend <- 6
  k_seasonal <- 6
  dat$nr_days <- 30
  f_null <- 'COUNT ~ s(EVENT_YEAR, k=k_trend, bs="ps")'
  if(!is.null(seasonal_spline_avg)){
    dat$avg_seasonal_val <- get_seasonal_spline_adj(seasonal_spline_avg,dat$EVENT_MONTH)
    offset <- 'offset(log(nr_days) + avg_seasonal_val)'
  }else{
    offset <- 'offset(log(nr_days))'
  }
  f_null <- paste(f_null,offset,sep='+')
  f_seasonal <- paste(f_null,'s(EVENT_MONTH, k=k_seasonal, bs="cp")',sep='+')
  mod_seasonal <- gam(as.formula(f_seasonal), family=quasipoisson(), data=dat, knots=list(EVENT_MONTH = c(a-0.5, b-0.5)), scale=-1)
  mod_null <- gam(as.formula(f_null), sp=mod_seasonal$sp['s(EVENT_YEAR)'],family=quasipoisson(), data=dat, scale=-1)
  return(list(seasonal=mod_seasonal, null=mod_null))
}

summarise_seasonality <- function(dat,seasonal_spline_avg=NULL){
  mod_list <- run_seasonality_gam(dat,seasonal_spline_avg)
  LR_test_stat <- anova(mod_list$null,mod_list$seasonal)
  log_pval <- pchisq(LR_test_stat$Deviance[2], df=LR_test_stat$Df[2], lower.tail=F, log.p=T)
  seasonal_component_monthly <- predict(mod_list$seasonal,
                                        newdata=filter(get_grid(dat,type='seasonal',seasonal_spline_avg),EVENT_MONTH %in% seq_len(12)),
                                        type='terms')[, 's(EVENT_MONTH)']
  seasonal_component_data <- predict(mod_list$seasonal, type='terms')[, 's(EVENT_MONTH)']
  val_peak <- exp(max(seasonal_component_monthly))
  val_trough <- exp(min(seasonal_component_monthly))
  month_peak <- which.max(seasonal_component_monthly)
  month_trough <- which.min(seasonal_component_monthly)
  ptr <- val_peak / val_trough
  log_relative_risk <- seasonal_component_monthly - min(seasonal_component_monthly)
  af <- group_by(dat,EVENT_MONTH) %>%
        summarise(count=sum(COUNT)) %>%
        ungroup() %>%
        mutate(p=count/sum(count),
               log_rr=log_relative_risk) %>%
        summarise(af=sum(p*(1-exp(-log_rr)))) %>%
        .$af
  tibble(pval=exp(log_pval),
         log10_pval=log_pval/log(10),
         month_peak=month_peak,
         month_trough=month_trough,
         val_peak=val_peak,
         val_trough=val_trough,
         ptr=ptr,
         af=af,
         var_fraction=max(0,1-var(mod_list$seasonal$residuals)/var(seasonal_component_data + mod_list$seasonal$residuals)),
         dispersion=mod_list$seasonal$scale)
}

extract_seasonal_spline <- function(dat){
  mod_list <- run_seasonality_gam(dat)
  months_vec <- seq(0.5,12.5,by=0.01)
  seasonal_pred <- predict(mod_list$seasonal,newdata=data.frame(nr_days=30,EVENT_YEAR=1998,EVENT_MONTH=months_vec),type='terms')[,'s(EVENT_MONTH)']
  return(tibble(month=months_vec,seasonal_val=seasonal_pred))
}

compare_gam_fits <- function(mod_list,monthly_counts,endpoint){
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
  ggplot(mutate(monthly_counts,month=1:n())) +
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
  scale_x_date(expand=c(0.005,0),date_breaks='1 year',date_labels='%Y') +
  scale_fill_manual(values=c('gold3','turquoise')) +
  xlab(paste0('First diagnosis of: ',endpoint)) +
  ylab('Count') +
  theme_classic() +
  theme(legend.position='none',
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14))
}

trend_plot <- function(dat,mod_list,seasonal_spline_avg=NULL){
  years_vec <- seq(1998,2019,by=0.1)
  trend_pred <- predict(mod_list$seasonal,newdata=get_grid(dat,type='trend',seasonal_spline_avg),type='terms',se.fit=T)
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

seasonality_plot <- function(dat,mod_list,a,b,seasonal_spline_avg=NULL){
  months_vec <- seq(a-0.5,b-0.5,by=0.1)
  seasonal_pred <- predict(mod_list$seasonal,newdata=get_grid(dat,type='seasonal',seasonal_spline_avg),type='terms',se.fit=T)
  seasonal_pred_dat <- tibble(month=months_vec,
                              est=seasonal_pred$fit[,'s(EVENT_MONTH)'],
                              lower=est - 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'],
                              upper=est + 1.96*seasonal_pred$se.fit[,'s(EVENT_MONTH)'])
  ggplot(seasonal_pred_dat) +
  geom_line(aes(x=month,
                y=est)) +
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper),alpha=0.3,fill='blue') +
  geom_hline(yintercept=0,linetype='dashed',color='red') +
  scale_x_continuous(breaks=seq(a,b)) +
  xlab('Month number') +
  ylab('Smooth term') +
  ggtitle('Seasonal component') +
  theme_bw() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14))
}
