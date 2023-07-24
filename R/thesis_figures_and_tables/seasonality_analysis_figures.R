library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(xtable)
library(mgcViz)
source('seasonality_dashboard/utils/seasonality_utils.R')
source('aux/tix_plotting.R')

#### Help functions ####
reorder_by_median <- function(dat,group_var,order='desc'){
  dat$group_var <- dat[[group_var]]
  median_dat <- group_by(dat,brief) %>%
                summarise(median_group=median(group_var))

  if(order=='desc'){
    median_dat <- arrange(median_dat,desc(median_group))
  }else{
    median_dat <- arrange(median_dat,median_group)
  }
  mutate(dat,brief=factor(brief,levels=median_dat$brief)) %>%
  select(-group_var)
}

split_brief <- function(s){
  ifelse(grepl(' and ',s),gsub(' and ',' and\n',s),gsub(' ','\n',s))
}


plot_seasonal_comp <- function(monthly_counts_FinRegistry,endpoint,endpoint_name){
  monthly_counts_endpoint <- filter(monthly_counts_FinRegistry,ENDPOINT==endpoint)
  unadj_model <- run_seasonality_gam(monthly_counts_endpoint,adjustment='',a=a,b=b)
  adj_model <- run_seasonality_gam(monthly_counts_endpoint,a=a,b=b,adjustment='binary',seasonal_spline_avg = avg_seasonal_spline)

  months_vec <- seq(a-0.5,b-0.5,by=0.1)
  seasonal_pred <- list(unadj=predict(unadj_model$seasonal,
                                      newdata=get_grid(type='seasonal'),
                                      type='terms',se.fit=T),
                        adj=predict(adj_model$seasonal,
                                    newdata=get_grid(type='seasonal',adjustment='binary'),
                                    type='terms',se.fit=T))

  seasonal_pred_dat <- tibble(model=c(rep('Unadjusted',length(months_vec)),rep('Adjusted',length(months_vec))),
                              month=rep(months_vec,2),
                              est=c(seasonal_pred$unadj$fit[,'s(EVENT_MONTH)'],seasonal_pred$adj$fit[,'s(EVENT_MONTH)']),
                              lower=est - 1.96*c(seasonal_pred$unadj$se.fit[,'s(EVENT_MONTH)'],seasonal_pred$adj$se.fit[,'s(EVENT_MONTH)']),
                              upper=est + 1.96*c(seasonal_pred$unadj$se.fit[,'s(EVENT_MONTH)'],seasonal_pred$adj$se.fit[,'s(EVENT_MONTH)']))
  seasonal_cmp_plot <- mutate(seasonal_pred_dat,model=factor(model,levels=c('Unadjusted','Adjusted'))) %>%
                        ggplot() +
                        geom_line(aes(x=month,
                                      y=est,col=model),linewidth=1) +
                        geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=model),alpha=0.3) +
                        geom_hline(yintercept=0,linetype='dashed',color='black') +
                        scale_x_continuous(breaks=seq(a,b-1,by=2)) +
                        scale_color_manual(values=c("#BC3C29FF","#0072B5FF"),name='Model type') +
                        scale_fill_manual(values=c("#BC3C29FF","#0072B5FF"),name='Model type') +
                        xlab('') +
                        ylab('Smooth term') +
                        theme_classic()+
                        theme(legend.position='none')


  years <- unique(monthly_counts_endpoint$EVENT_YEAR)
  season_dat <- bind_rows(list(summer=tibble(start=ym(paste0(years,'-',5)),
                                             end=ym(paste0(years,'-',9))),
                               winter=tibble(start=ym(paste0(years,'-',10)),
                                             end=ym(paste0(years+1,'-',4)))),.id='season')
  season_dat <- bind_rows(season_dat,tibble(season='winter',
                                            start=ym(paste0(years[1],'-',1)),
                                            end=ym(paste0(years[1],'-',4))))
  monthly_counts_endpoint$null_pred <-predict(unadj_model$null)
  monthly_counts_endpoint$seasonal_pred <-predict(unadj_model$seasonal)
  monthly_counts_endpoint <- mutate(monthly_counts_endpoint,month=1:n())
  fit_plot <- ggplot(monthly_counts_endpoint) +
              geom_point(aes(x=EVENT_DATE,
                             y=COUNT),alpha=0.5) +
              geom_line(aes(x=EVENT_DATE,
                            y=exp(seasonal_pred))) +
              geom_line(aes(x=EVENT_DATE,
                            y=exp(null_pred)),
                        col='red') +
              geom_rect(data=season_dat,aes(xmin=start-50,
                                            xmax=end,
                                            ymin=0.97*min(monthly_counts_endpoint$COUNT),
                                            ymax=1.03*max(monthly_counts_endpoint$COUNT),
                                            fill=season),alpha=0.1) +
              scale_y_continuous(expand=c(0,0)) +
              scale_x_date(expand=c(0.005,0),date_breaks='3 years',date_labels='%Y') +
              scale_fill_manual(values=c('gold3','turquoise')) +
              ggtitle(endpoint_name) +
              xlab('') +
              ylab('Count') +
              theme_classic() +
              theme(legend.position='none')
  return(list(fit_plot=fit_plot,seasonal_cmp_plot=seasonal_cmp_plot))
}
###### Read necessary data and define endpoint categories #######
a <- 1;b <- 13
GWAS_info_dat <- read_tsv('data/finngen_R10_finngen_R10_analysis_data_finngen_R10_pheno_n.tsv')

monthly_counts_FinRegistry_org  <- read_tsv('data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
                               filter_monthly_counts(.,GWAS_info_dat=GWAS_info_dat,endpoint_mapping_path='data/finngen_R10_endpoint_core_noncore_1.0.xlsx')
endpoint_categories <- read_excel('data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=2)
endpoint_table <- read_excel('data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=1) %>%
                  filter(CORE_ENDPOINTS=='yes') %>%
                    select(NAME,TAGS,LONGNAME) %>%
                    separate_rows(TAGS,sep=',') %>%
                    left_join(endpoint_categories,by=c('TAGS'='code')) %>%
                    select(NAME,TAGS,LONGNAME,brief)

monthly_counts_FinRegistry  <- left_join(monthly_counts_FinRegistry_org,endpoint_table,by=c('ENDPOINT'='NAME'),suffix=c('_prefix','_tag')) %>%
                               mutate(brief_prefix=ifelse(is.na(brief_prefix),brief_tag,brief_prefix),
                                      brief_tag=ifelse(is.na(brief_tag),LONGNAME,brief_tag),
                                      brief=ifelse(brief_prefix!=brief_tag,brief_prefix,brief_tag))

seasonal_summary <- read_tsv('data/FINREGISTRY_seasonal_summary.txt')
seasonal_summary_binary_adj <- read_tsv('data/FINREGISTRY_seasonal_summary_binary_adj.txt')
seasonal_summary_mean_covariate_adj <- read_tsv('data/FINREGISTRY_seasonal_summary_mean_covariate_adj.txt')
seasonal_summary_mean_offset_adj <- read_tsv('data/FINREGISTRY_seasonal_summary_mean_offset_adj.txt')
seasonal_splines <- read_tsv('data/FINREGISTRY_seasonal_splines.txt')
seasonal_splines_adj <- read_tsv('data/FINREGISTRY_seasonal_splines_adj_binary.txt')


endpoint_mapping <- distinct(monthly_counts_FinRegistry,ENDPOINT,brief,prefix)
pval_thres <- 0.05/nrow(seasonal_summary)

#####  Supplementary table - Endpoint categories ######
endpoint_categories <- inner_join(seasonal_summary,endpoint_mapping,by='ENDPOINT') %>%
                       inner_join(seasonal_summary_binary_adj,by='ENDPOINT',suffix=c('','_adj')) %>%
                       group_by(brief) %>%
                       summarise(n=n(),
                                 n_signif=sum(smooth_pval<pval_thres),
                                 n_signif_adj=sum(smooth_pval_adj<pval_thres)) %>%
                       filter(n_signif>20) %>%
                       arrange(desc(n))
print(xtable(select(endpoint_categories,`Endpoint category`=brief,`Nr of endpoints`=n,`Nr of significant endpoints`=n_signif)),include.rownames=FALSE)



#### Seasonality overview figure #####
### Avg spline plot
avg_seasonal_spline <- group_by(seasonal_splines,month) %>%
                           summarise(avg_seasonal_val=mean(seasonal_val),
                                     median_seasonal_val=median(seasonal_val))
set.seed(10)
random_endpoints <- seasonal_summary$ENDPOINT[sample(1:nrow(seasonal_summary),size = 20)]
avg_spline_plot <- filter(seasonal_splines,ENDPOINT %in% random_endpoints) %>%
                    ggplot() +
                    geom_line(data=,aes(month,seasonal_val,group=ENDPOINT),col='black',alpha=0.2) +
                    geom_line(data=avg_seasonal_spline,aes(month,avg_seasonal_val),col='red',linewidth=1) +
                    scale_x_continuous(expand=c(0,0),breaks=seq_len(12)) +
                    scale_y_continuous(limits=c(-0.25,0.25),expand=c(0,0)) +
                    theme_bw() +
                    xlab('Month number') +
                    ylab(expression(s[s](x)))




#Delta PTR plots
PTR_percentage_dat <- inner_join(seasonal_summary,seasonal_summary_binary_adj,by='ENDPOINT',suffix=c('','_adj')) %>%
                        #filter(ptr_est>1.1) %>%
                        filter(smooth_pval<pval_thres) %>%
                        mutate(delta_PTR=-100*((ptr_est-ptr_est_adj)/(ptr_est-1))) %>%
                        select(ENDPOINT,delta_PTR) %>%
                        inner_join(endpoint_mapping,by='ENDPOINT') %>%
                        filter(!is.na(brief)) %>%
                        group_by(brief) %>%
                        filter(n()>20) %>%
                        ungroup() %>%
                        mutate(brief=ifelse(nchar(brief)>20,split_brief(brief),brief)) %>%
                        #filter(abs(delta_PTR)<1) %>%
                        reorder_by_median('delta_PTR',order='desc')

PTR_percentage_dat %>% group_by(brief) %>% summarise(median(delta_PTR))

PTR_percentage_plot <- ggplot(PTR_percentage_dat,aes(brief,y=delta_PTR)) +
                        geom_boxplot(alpha=0.5,fill="#0072B5FF") +
                        geom_hline(yintercept=0,color='red') +
                        #scale_y_continuous(limits=c(0,NA),expand=c(0.01,0)) +
                        xlab('') +
                        ylab('Delta PTR \\%') +
                        theme_classic() +
                        theme(axis.text.x = element_text(size=8,angle=-45,hjust=0),
                              plot.margin = margin(0,3,0,0, "cm"))

#season peak
season_peak_dat <- filter(seasonal_summary,smooth_pval<pval_thres) %>%
                    inner_join(endpoint_mapping,by='ENDPOINT') %>%
                    mutate(season_peak=case_when(month_peak<5 ~ 'winter',
                                                 month_peak<9 ~ 'summer',
                                                 TRUE ~ 'winter')) %>%
                    filter(!is.na(brief)) %>%
                    group_by(brief) %>%
                    filter(n()>20) %>%
                    ungroup() %>%
                    count(brief,season_peak,name = 'count') %>%
                    left_join(expand_grid(brief=unique(.$brief),season_peak=c('summer','winter')),.,by=c('brief','season_peak')) %>%
                    mutate(count=ifelse(is.na(count),0,count)) %>%
                    group_by(brief) %>%
                    mutate(pct=100*count/sum(count)) %>%
                    ungroup() %>%
                    mutate(brief=ifelse(nchar(brief)>20,split_brief(brief),brief))
filter(season_peak_dat,grepl('Infections|Pregnancy',brief))


season_peak_plot <- ggplot(season_peak_dat,aes(brief,y=count,fill=season_peak)) +
                    geom_col(position='dodge',alpha=0.5) +
                    scale_y_continuous(expand=c(0,0)) +
                    scale_fill_manual(values=c("#E18727FF","#0072B5FF"),name='Season peak') +
                    xlab('') +
                    ylab('Nr of endpoints') +
                    theme_classic() +
                    theme(axis.text.x = element_text(size=8,angle=-45,hjust=0),
                          legend.position=c(0.8,0.8),
                          plot.margin = margin(0,3,0,0, "cm"))

#month peak
month_peak_dat <- filter(seasonal_summary,(smooth_pval)<pval_thres) %>%
                  mutate(month_peak=ifelse(month_peak<1,floor(month_peak+12),floor(month_peak))) %>%
                  count(month_peak) %>%
                  mutate(pct=100*n/sum(n))

filter(month_peak_dat,month_peak>=4 & month_peak<=8) %>%
summarise(sum(pct))

filter(month_peak_dat,month_peak %in% c(2,3,10)) %>%
summarise(sum(pct))

month_peak_plot <- ggplot(month_peak_dat,aes(month_peak,y=pct)) +
                    geom_col(position='dodge',fill="#BC3C29FF",alpha=0.5,width=0.8) +
                    scale_x_continuous(breaks=seq(1,12),expand=c(0,0.2)) +
                    scale_y_continuous(expand=c(0,0)) +
                    xlab('Peak month number') +
                    ylab('\\% of endpoints') +
                    theme_classic() +
                    theme(plot.margin = margin(0,0,0,0, "cm"),
                          axis.title.x = element_text(vjust = 1))


seasonal_overview_plot <- ((  + ggtitle('A'))  + (month_peak_plot + ggtitle('B'))) /
                          ((season_peak_plot + ggtitle('C'))) /
                          ((PTR_percentage_plot + ggtitle('D')))
save_tikz_plot(seasonal_overview_plot,width=7.2,height=7.5,filename = 'tex/seasonal_overview_plot.tex')

##### Seasonality of four endpoint with reported seasonality in the literature #####


seasonal_autoimmune=plot_seasonal_comp(monthly_counts_FinRegistry_org,endpoint='AUTOIMMUNE',endpoint_name='Autoimmune disease')
#seasonal_cmp_acne <- plot_seasonal_comp(monthly_counts_FinRegistry,endpoint='L12_ACNE',endpoint_name='Acne')
seasonal_influenza <- plot_seasonal_comp(monthly_counts_FinRegistry_org,endpoint='J10_INFLUENZA',endpoint_name='Influenza')
seasonal_depression <- plot_seasonal_comp(monthly_counts_FinRegistry_org,endpoint='F5_DEPRESSIO',endpoint_name='Major depression')
seasonal_CVD <- plot_seasonal_comp(monthly_counts_FinRegistry_org,endpoint='FG_CVD',endpoint_name='Cardiovascular disease')

seasonal_plots_disease <- seasonal_autoimmune$fit_plot + seasonal_autoimmune$seasonal_cmp_plot +
                          seasonal_influenza$fit_plot + seasonal_influenza$seasonal_cmp_plot +
                          (seasonal_depression$fit_plot) + (seasonal_depression$seasonal_cmp_plot) + #theme(legend.position=c(0.6,0.8))) +
                          (seasonal_CVD$fit_plot + xlab('First diagnosis date')) + (seasonal_CVD$seasonal_cmp_plot + xlab('Month number')) +
                          plot_layout(ncol=2,nrow=4,byrow=T,widths = c(16, 4))
save_tikz_plot(seasonal_plots_disease,width=6.5,height=7.5,filename = 'tex/seasonality_plot_diseases.tex')

##### Summary table of the four endpoints with reported seasonality in the literature #####
### Seasonal summary table
endpoint_name_mapping <- read_excel('R/thesis_figures_and_tables/endpoint_name_mapping.xlsx')
endpoints <- c('AUTOIMMUNE','F5_DEPRESSIO','FG_CVD','J10_INFLUENZA')

adj_vs_undaj <- bind_rows(
  list(adj=seasonal_summary_binary_adj %>% filter(ENDPOINT %in% endpoints),
       unadj=seasonal_summary %>% filter(ENDPOINT %in% endpoints)),.id = 'type') %>%
  select(ENDPOINT,type,ptr_est,ptr_lower,ptr_upper,month_peak,month_peak_lower,month_peak_upper,month_trough,month_trough_lower,month_trough_upper) %>%
  mutate(type=factor(type,levels=c('unadj','adj'),labels=c('Unadjusted','Adjusted'))) %>%
arrange(ENDPOINT,type) %>%
#pivot_wider(id_cols='ENDPOINT',names_from = 'type',values_from=c('ptr_est','ptr_lower','ptr_upper','month_peak','month_trough')) %>%
mutate(ptr_est=paste0(formatC(ptr_est,format='f',digits=2),' (',formatC(ptr_lower,format='f',digits=2),'-',formatC(ptr_upper,format='f',digits=2),')'),
       month_peak=paste0(formatC(month_peak,format='f',digits=2),' (',formatC(month_peak_lower,format='f',digits=2),'-',formatC(month_peak_upper,format='f',digits=2),')'),
       month_trough=paste0(formatC(month_trough,format='f',digits=2),' (',formatC(month_trough_lower,format='f',digits=2),'-',formatC(month_trough_upper,format='f',digits=2),')')) %>%
inner_join(endpoint_name_mapping,by=c('ENDPOINT'='Endpoint ID')) %>%
mutate(Endpoint=ifelse(type=='Adjusted','',Endpoint)) %>%
select(Endpoint,type,ptr_est,month_peak,month_trough)


print(xtable(select(adj_vs_undaj,Endpoint,Model=type,,PTR=ptr_est,`Month peak`=month_peak,`Month trough`=month_trough)),include.rownames=FALSE)

###### Seasonality figure for heart failure ######
seasonal_heartfail <- plot_seasonal_comp(monthly_counts_FinRegistry,endpoint='I9_HEARTFAIL_ALLCAUSE',endpoint_name='Heart failure')

#Heart failure
seasonal_plot_heartfail <-(seasonal_heartfail$fit_plot + xlab('First diagnosis date')) +
                          (seasonal_heartfail$seasonal_cmp_plot + xlab('Month number')) +# + theme(legend.position=c(0.8,0.8))) +
                          plot_layout(ncol=2,nrow=1,byrow=T,widths = c(16, 4))

save_tikz_plot(seasonal_plot_heartfail,width=6.5,height=3,filename = 'tex/seasonality_plot_heartfail.tex')

####### Example of model diagnostics for CVD #####
#Model diagnostics
mod_CVD=run_seasonality_gam(monthly_counts_FinRegistry_org %>% filter(ENDPOINT=='FG_CVD'),a=1,b=13)
check(getViz(mod_CVD$seasonal))

##### Visualization of odispersion parameter #####
#Over-dispersion
dispersion_plot <- ggplot(seasonal_summary,aes(dispersion)) +
                    geom_histogram(bins=50,col='black',fill='turquoise',alpha=0.5) +
                    scale_x_log10(expand=c(0,0)) +
                    scale_y_continuous(expand=c(0.001,0)) +
                    xlab('Dispersion (log10 scale)') +
                    ylab('Count') +
                    theme_classic()
save_tikz_plot(dispersion_plot,width=6,height=4,filename = 'tex/dispersion_hist.tex')

count(seasonal_summary,dispersion>1)




### Deviance explained for different adjustment models
dev_dat <- inner_join(select(seasonal_summary_mean_offset_adj,ENDPOINT,deviance,deviance_expl),
           select(seasonal_summary_mean_covariate_adj,ENDPOINT,deviance,deviance_expl),by='ENDPOINT',suffix=c('_mean_offset','_mean_covariate')) %>%
           inner_join(select(seasonal_summary_binary_adj,ENDPOINT,deviance,deviance_expl),by='ENDPOINT')


dev_diff_dat <- mutate(dev_dat,mean_offset_diff=deviance_expl-deviance_expl_mean_offset,
                               mean_covariate_diff=deviance_expl-deviance_expl_mean_covariate) %>%
                select(ENDPOINT,mean_offset_diff,mean_covariate_diff) %>%
                pivot_longer(cols=c('mean_offset_diff','mean_covariate_diff')) %>%
                mutate(name=factor(name,levels=c('mean_offset_diff','mean_covariate_diff'),
                                        labels=c('Mean seasonality curve offset','Mean seasonality curve covariate')))

adj_deviance_diff_plot <- ggplot(dev_diff_dat,aes(name,value)) +
                          geom_boxplot(width=0.4,fill='turquoise',alpha=0.7) +
                          geom_violin(width=0.4,fill='turquoise',alpha=0.3) +
                          geom_hline(yintercept=0) +
                          xlab('') +
                          ylab('Deviance explained difference') +
                          theme_classic()
save_tikz_plot(adj_deviance_diff_plot,width=6,height=4,filename = 'tex/adj_deviance_diff_plot.tex')




# dev_diff_seasonality_dat <- mutate(dev_dat,mean_offset_diff=deviance_expl_seasonality-deviance_expl_seasonality_mean_offset,
#                                    mean_covariate_diff=deviance_expl_seasonality-deviance_expl_seasonality_mean_covariate) %>%
#                             select(ENDPOINT,mean_offset_diff,mean_covariate_diff) %>%
#                             pivot_longer(cols=c('mean_offset_diff','mean_covariate_diff')) %>%
#                             mutate(name=factor(name,levels=c('mean_offset_diff','mean_covariate_diff'),
#                                                labels=c('Mean seasonality curve offset','Mean seasonality curve covariate')))
#
#
# adj_deviance_diff_seasonality_plot <- ggplot(dev_diff_seasonality_dat,aes(name,value)) +
#                                       geom_boxplot(width=0.4,fill='turquoise',alpha=0.7) +
#                                       geom_violin(width=0.4,fill='turquoise',alpha=0.3) +
#                                       geom_hline(yintercept=0) +
#                                       xlab('') +
#                                       ylab('Deviance explained difference') +
#                                       theme_classic() +
#                                       theme_custom



#
#
#   bind_rows(list(
#     unadj=filter(seasonal_summary,pval<pval_thres),
#     adj=filter(seasonal_summary_binary_adj,pval<pval_thres)),.id='mod_type') %>%
#   mutate(mod_type=factor(mod_type,levels=c('unadj','adj'))) %>%
#   inner_join(endpoint_mapping,by='ENDPOINT') %>%
#   count(prefix,descr,brief,mod_type) %>%
#   ungroup() %>%
#   reorder_by_count(filter_var='unadj') %>%
#   ggplot() +
#   geom_col(aes(x=brief,y=n,fill=mod_type),position='dodge') +
#   scale_y_continuous(expand=c(0,0)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=-30,hjust=0),
#         legend.position=c(0.8,0.8),
#         plot.margin = margin(0.5,6,0,0.5, "cm"))
#
#
# filter(seasonal_summary_binary_adj,val_peak<10) %>%
#   ggplot(aes(log(val_peak))) +
#   geom_histogram(bins=100,fill='red',alpha=0.3,col='black') +
#   scale_x_continuous(expand=c(0,0.01)) +
#   scale_y_continuous(expand=c(0,0.5)) +
#   xlab('Log of seasonal component peak') +
#   ylab('Count') +
#   theme_classic() +
#   theme_custom
#
# filter(seasonal_summary,val_peak<10) %>%
#   ggplot(aes(log(val_trough))) +
#   geom_histogram(bins=100,fill='blue',alpha=0.3,col='black') +
#   scale_x_continuous(expand=c(0,0.01)) +
#   scale_y_continuous(expand=c(0,0.5)) +
#   xlab('Log of seasonal component trough') +
#   ylab('Count') +
#   theme_classic() +
#   theme_custom

### Hospital visits plots
#hospital_visits_weights <- read_tsv('seasonality_dashboard/data/hospital_visits_weights.txt')

# hospital_visits_plot <- ggplot(hospital_visits_weights,aes(EVENT_MONTH,100*COUNT_PROP)) +
#   geom_col(fill='#BC3C29FF',alpha=0.5) +
#   scale_y_continuous(expand=c(0,0),limits=c(0,11),breaks=seq(0,10,by=2)) +
#   scale_x_continuous(expand=c(0.01,0),breaks = seq(1,12)) +
#   xlab('Month number') +
#   ylab('\\% of hospital visits') +
#   theme_classic() +
#   theme(plot.margin = margin(0,0,0,0, "cm"))

# plot_cmp_seasonal_metric <- function(seasonal_summary,seasonal_summary_adj,endpoint_mapping,var,transformation=identity){
#   seasonal_summary$var_transformed <- transformation(seasonal_summary[[var]])
#   seasonal_summary <- filter(seasonal_summary,val_peak<5) %>%
#     inner_join(endpoint_mapping,by='ENDPOINT') %>%
#     filter(!is.na(prefix)) %>%
#     group_by(prefix) %>%
#     filter(n()>20) %>%
#     ungroup() %>%
#     reorder_by_median('var_transformed')
#
#   seasonal_summary_adj$var_transformed <- transformation(seasonal_summary_adj[[var]])
#   seasonal_summary_adj <- filter(seasonal_summary_adj,val_peak<5) %>%
#     inner_join(endpoint_mapping,by='ENDPOINT') %>%
#     filter(!is.na(prefix)) %>%
#     group_by(prefix) %>%
#     filter(n()>20) %>%
#     ungroup() %>%
#     mutate(brief=factor(brief,levels=levels(seasonal_summary$brief)),
#            prefix=factor(prefix,levels=levels(seasonal_summary$prefix)))
#   bind_rows(list(Unadjusted=seasonal_summary,
#                  Adjusted=seasonal_summary_adj),.id='model_type') %>%
#     mutate(model_type=factor(model_type,levels=c('Unadjusted','Adjusted')))  %>%
#     ggplot(aes(brief,y=var_transformed,fill=model_type)) +
#     geom_boxplot(alpha=0.6) +
#     scale_fill_manual(values=c("#BC3C29FF","#0072B5FF"),name='Model type') +
#     xlab('') +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle=-30,hjust=0),
#           plot.margin = margin(0,3,0,0, "cm"))
# }

# val_peak_plot <- plot_cmp_seasonal_metric(seasonal_summary,seasonal_summary_binary_adj,endpoint_mapping,var='val_peak',transformation=log) +
#                  scale_y_continuous(limits=c(0,0.5),breaks=seq(0,0.5,by=0.1),expand=c(0,0.02)) +
#                  ylab('Log of peak value') +
#                  theme(legend.position=c(0.8,0.8))

# var_fraction_plot <- plot_cmp_seasonal_metric(seasonal_summary,seasonal_summary_binary_adj,endpoint_mapping,var='var_fraction') +
#                      scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.1),expand=c(0,0.02)) +
#                      ylab('Variance explained by season') +
#                      theme(legend.position=c(0.8,0.8))

# reorder_by_count <- function(dat,filter_var){
#   count_dat <- filter(dat,mod_type==filter_var) %>%
#     arrange(desc(n))
#
#   out_dat <- mutate(dat,prefix=factor(prefix,levels=count_dat$prefix),
#                     brief=factor(brief,levels=count_dat$brief))
#   return(out_dat)
# }

# plot_seasonal_metric <- function(seasonal_summary,var,transformation=identity){
#   seasonal_summary$var_transformed <- transformation(seasonal_summary[[var]])
#   filter(seasonal_summary,val_peak<5) %>%
#     inner_join(endpoint_mapping,by='ENDPOINT') %>%
#     group_by(prefix) %>%
#     filter(n()>20) %>%
#     ungroup() %>%
#     reorder_by_median('var_transformed')%>%
#     ggplot(aes(brief,y=var_transformed)) +
#     geom_boxplot(fill='red',alpha=0.3) +
#     theme_classic() +
#     xlab('') +
#     theme(axis.text.x = element_text(angle=-30,hjust=0),
#           plot.margin = margin(0.5,6,0,0.5, "cm"))
# }



