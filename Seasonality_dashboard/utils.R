run_seasonality_gam <- function(dat){
  k_trend <- 11
  k_seasonal <- 6
  dat$nr_days <- 30
  f_null <- 'COUNT ~ offset(log(nr_days)) + s(EVENT_YEAR, k=k_trend, bs="ps")'
  f_seasonal <- paste(f_null,'s(EVENT_MONTH, k=k_seasonal, bs="cp")',sep='+')
  mod_null <- gam(as.formula(f_null), family=quasipoisson(), data=dat, scale=-1)
  mod_seasonal <- gam(as.formula(f_seasonal), family=quasipoisson(), data=dat, scale=-1)
  return(list(seasonal=mod_seasonal, null=mod_null))
}

summarise_seasonality <- function(dat){
  mod_list <- run_seasonality_gam(dat)
  LR_test_stat <- anova(mod_list$null,mod_list$seasonal)
  log_pval <- pchisq(LR_test_stat$Deviance[2], df=LR_test_stat$Df[2], lower.tail=F, log.p=T)
  seasonal_component_monthly <- predict(mod_list$seasonal,
                                        newdata=data.frame(nr_days=30, EVENT_YEAR=1998, EVENT_MONTH=seq_len(12)),
                                        type='terms')[, 2]
  seasonal_component_data <- predict(mod_list$seasonal, type='terms')[, 2]
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

compare_gam_fits <- function(mod_list,monthly_counts,endpoint){
  years <- unique(monthly_counts$EVENT_YEAR)
  season_dat <- bind_rows(list(summer=tibble(start=ym(paste0(years,'-',4)),
                                             end=ym(paste0(years,'-',8))),
                               winter=tibble(start=ym(paste0(years,'-',8)),
                                             end=ym(paste0(years+1,'-',4)))),.id='season')
  season_dat <- bind_rows(season_dat,tibble(season='winter',
                                            start=ym(paste0(years[1],'-',1)),
                                            end=ym(paste0(years[1],'-',4))))
  ggplot(mutate(monthly_counts,month=1:n())) +
  geom_point(aes(x=EVENT_DATE,
                 y=COUNT)) +
  geom_line(aes(x=EVENT_DATE,
                y=mod_list$seasonal$fitted.values)) +
    geom_line(aes(x=EVENT_DATE,
                  y=mod_list$null$fitted.values),
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

trend_plot <- function(mod_list){
  years_vec <- seq(1998,2019,by=0.1)
  trend_vec <- predict(mod_list$seasonal,newdata=data.frame(nr_days=30,EVENT_YEAR=years_vec,EVENT_MONTH=1),type='terms')[,1]
  tibble(year=years_vec,trend=trend_vec) %>%
  ggplot(aes(x=year,
             y=trend)) +
  geom_line() +
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

seasonality_plot <- function(mod_list){
  months_vec <- seq(1,12,by=0.1)
  seasonal_vec <- predict(mod_list$seasonal,newdata=data.frame(nr_days=30,EVENT_YEAR=1998,EVENT_MONTH=months_vec),type='terms')[,2]
  tibble(month=months_vec,seasonal=seasonal_vec) %>%
  ggplot(aes(x=months_vec,
             y=seasonal)) +
  geom_line() +
  scale_x_continuous(breaks=seq_len(12)) +
  xlab('Month number') +
  ylab('Smooth term') +
  ggtitle('Seasonal component') +
  theme_bw() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=14))
}


filter_monthly_counts <- function(monthly_counts,heritability_dat){
  monthly_counts_filtered <- filter(monthly_counts,EVENT_YEAR>=1998 & EVENT_YEAR<=2019) %>%
                             inner_join(select(heritability_dat,ENDPOINT=PHENO,H2),by='ENDPOINT') %>%
                             filter(H2>0.01)
  month_year_endpoint_template <- expand_grid(ENDPOINT=unique(monthly_counts_filtered$ENDPOINT),
                                              EVENT_YEAR=seq(min(monthly_counts_filtered$EVENT_YEAR),
                                                             max(monthly_counts_filtered$EVENT_YEAR)),
                                              EVENT_MONTH=seq_len(12))
  monthly_counts_filtered <- left_join(month_year_endpoint_template,
                                       monthly_counts_filtered,
                                       by=c('ENDPOINT','EVENT_YEAR','EVENT_MONTH')) %>%
                              mutate(COUNT=ifelse(is.na(COUNT),0L,COUNT),
                                     EVENT_DATE=ym(paste0(EVENT_YEAR,'-',EVENT_MONTH))) %>%
                              arrange(ENDPOINT,EVENT_YEAR,EVENT_MONTH)
  return(monthly_counts_filtered)
}


# t = Sys.time()
# seasonal_summary_FinRegistry = group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
#                                group_modify(~summarise_seasonality(.x)) %>%
#                                arrange(log10_pval)
# Sys.time()-t
# write.table(seasonal_summary_FinRegistry,file='data/FINREGISTRY_seasonal_summary.txt',sep='\t',row.names=F,quote=F)

