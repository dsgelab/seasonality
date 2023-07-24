library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(patchwork)
library(xtable)
library(gseasonality)
source('seasonality_paper/R/utils.R')
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

###### Read necessary data and define endpoint categories #######
seasonality_gam_FinRegistry_unadj <- readRDS('seasonality_paper/results/seasonality_gam_FinRegistry_unadj_1999-2019.rds')
seasonality_gam_FinRegistry_adj <- readRDS('seasonality_paper/results/seasonality_gam_FinRegistry_adj_1999-2019.rds')

seasonality_gam_summary_table <- lapply(names(seasonality_gam_FinRegistry_unadj),function(e){
                                    tibble(ENDPOINT=e,
                                           ptr_est_unadj=seasonality_gam_FinRegistry_unadj[[e]]$summary$ptr$estimate,
                                           ptr_est_adj=seasonality_gam_FinRegistry_adj[[e]]$summary$ptr$estimate,
                                           seasonality_pval_unadj=seasonality_gam_FinRegistry_unadj[[e]]$summary$pval,
                                           seasonality_pval_adj=seasonality_gam_FinRegistry_adj[[e]]$summary$pval,
                                           month_peak_unadj=seasonality_gam_FinRegistry_unadj[[e]]$summary$peak_trough['month','peak'],
                                           month_peak_adj=seasonality_gam_FinRegistry_adj[[e]]$summary$peak_trough['month','peak'],
                                           month_trough_unadj=seasonality_gam_FinRegistry_unadj[[e]]$summary$peak_trough['month','trough'],
                                           month_trough_adj=seasonality_gam_FinRegistry_adj[[e]]$summary$peak_trough['month','trough'])
                                  }) %>% bind_rows()

seasonal_spline_avg_unadj <- lapply(names(seasonality_gam_FinRegistry_unadj),function(e) mutate(ENDPOINT=e,seasonality_gam_FinRegistry_unadj[[e]]$seasonality_term)) %>%
                             bind_rows() %>%
                             group_by(month) %>%
                             summarise(avg_seasonal_val=mean(est))

seasonal_spline_avg_adj <- lapply(names(seasonality_gam_FinRegistry_adj),function(e) mutate(ENDPOINT=e,seasonality_gam_FinRegistry_adj[[e]]$seasonality_term)) %>%
                            bind_rows() %>%
                            group_by(month) %>%
                            summarise(avg_seasonal_val=mean(est))

endpoint_categories <- read_excel('seasonality_paper/data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=2)
endpoint_mapping_manual <- read_excel('seasonality_paper/data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=1) %>%
                            filter(CORE_ENDPOINTS=='yes') %>%
                            select(NAME,TAGS,LONGNAME) %>%
                            separate_rows(TAGS,sep=',') %>%
                            left_join(endpoint_categories,by=c('TAGS'='code')) %>%
                            select(NAME,TAGS,LONGNAME,brief)
endpoint_mapping_prefix <- get_endpoint_mapping(names(seasonality_gam_FinRegistry_unadj),endpoint_categories)
endpoint_mapping  <- left_join(endpoint_mapping_prefix,endpoint_mapping_manual,by=c('ENDPOINT'='NAME'),suffix=c('_prefix','_manual')) %>%
                     mutate(brief_prefix=ifelse(is.na(brief_prefix),brief_manual,brief_prefix),
                            brief_manual=ifelse(is.na(brief_manual),LONGNAME,brief_manual),
                            brief=ifelse(is.na(brief_manual),brief_prefix,brief_manual),
                            brief=ifelse(brief_manual!=brief_prefix,brief_prefix,brief_manual)) %>%
                    distinct(ENDPOINT,brief,prefix) %>%
                    group_by(brief) %>%
                    filter(n()>20) %>%
                    ungroup() %>%
                    filter(!is.na(brief))

pval_thres <- 0.05/length(seasonality_gam_FinRegistry_unadj)

#####  Supplementary table - Endpoint categories ######
endpoint_categories <- inner_join(seasonality_gam_summary_table,endpoint_mapping,by='ENDPOINT') %>%
                      group_by(brief) %>%
                      summarise(n=n(),
                                n_signif_unadj=sum(seasonality_pval_unadj<pval_thres),
                                n_signif_adj=sum(seasonality_pval_adj<pval_thres)) %>%
                      filter(n_signif_unadj>20) %>%
                      arrange(desc(n))
print(xtable(select(endpoint_categories,`Endpoint category`=brief,`Nr of endpoints`=n,`Nr of significant endpoints`=n_signif_unadj)),include.rownames=FALSE)



#### Figure 2 - Epidemiological overview #####
### Avg spline plot
set.seed(10)
random_endpoints_seasonality_term <- lapply(sample(1:length(seasonality_gam_FinRegistry_unadj),size=20,replace=F),function(i){
                                        mutate(seasonality_gam_FinRegistry_unadj[[i]]$seasonality_term,
                                               ENDPOINT=names(seasonality_gam_FinRegistry_unadj)[i])
                                      }) %>% bind_rows()
avg_spline_plot <- ggplot() +
                    geom_line(data=random_endpoints_seasonality_term,aes(month,exp(est),group=ENDPOINT),col='black',alpha=0.2) +
                    geom_line(data=seasonal_spline_avg_unadj,aes(month,exp(avg_seasonal_val)),col='red',linewidth=1) +
                    scale_x_continuous(expand=c(0,0),breaks=seq_len(12)) +
                    #scale_y_continuous(limits=c(-0.25,0.25),expand=c(0,0)) +
                    theme_bw() +
                    xlab('Month number') +
                    ylab('Seasonality term')


peak_trough_hist <- seasonality_gam_summary_table %>%
                    filter(seasonality_pval_unadj<pval_thres) %>%
                    mutate(month_peak_unadj_rounded=round(ifelse(month_peak_unadj<1,month_peak_unadj+12,month_peak_unadj),1),
                           month_trough_unadj_rounded=round(ifelse(month_trough_unadj<1,month_trough_unadj+12,month_trough_unadj),1)) %>%
                    pivot_longer(cols=c('month_peak_unadj_rounded','month_trough_unadj_rounded')) %>%
                    group_by(name,value) %>%
                    summarise(n=n()) %>%
                    mutate(perc=100*n/sum(n)) %>%
                    ungroup() %>%
                    mutate(perc=ifelse(grepl('trough',name),-perc,perc),
                           name=factor(name,levels=c('month_peak_unadj_rounded','month_trough_unadj_rounded'),labels=c('Peak','Trough'))) %>%
                    ggplot() +
                    geom_col(aes(x=value,y=perc,fill=name),alpha=1) +
                    geom_hline(yintercept=0) +
                    scale_x_continuous(breaks=1:12) +
                    scale_y_continuous(breaks=seq(-20,10,by=10),labels = c(seq(20,0,by=-10),10),limits=c(-25,10)) +
                    scale_fill_manual(values=c('#BC3C29FF','#0072B5FF'),name='') +
                    xlab('Month number') +
                    ylab('% of endpoints') +
                    theme_minimal() +
                    theme(legend.position = c(0.8,0.4))

# peak_trough_hist <- seasonality_gam_summary_table %>%
#                     filter(seasonality_pval_unadj<pval_thres) %>%
#                     mutate(month_peak_unadj_int=floor(ifelse(month_peak_unadj<1,month_peak_unadj+12,month_peak_unadj)),
#                            month_trough_unadj_int=floor(ifelse(month_trough_unadj<1,month_trough_unadj+12,month_trough_unadj))) %>%
#                     pivot_longer(cols=c('month_peak_unadj_int','month_trough_unadj_int')) %>%
#                     group_by(name,value) %>%
#                     summarise(n=n()) %>%
#                     mutate(perc=100*n/sum(n)) %>%
#                     ungroup() %>%
#                     mutate(perc=ifelse(grepl('trough',name),-perc,perc),
#                            name=factor(name,levels=c('month_peak_unadj_int','month_trough_unadj_int'),labels=c('Peak','Trough'))) %>%
#                     ggplot() +
#                     geom_col(aes(x=value,y=perc,fill=name),col='black',alpha=0.7) +
#                     geom_hline(yintercept=0) +
#                     scale_x_continuous(breaks=1:12) +
#                     scale_y_continuous(breaks=seq(-60,40,by=20),labels = c(seq(60,0,by=-20),20,40),limits=c(-75,35)) +
#                     scale_fill_manual(values=c('#BC3C29FF','#0072B5FF'),name='') +
#                     xlab('Month number') +
#                     ylab('\\% of endpoints') +
#                     theme_minimal() +
#                     theme(legend.position = c(0.8,0.3))

#season peak
# season_peak_dat <- filter(seasonality_gam_summary_table,seasonality_pval_adj<pval_thres) %>%
#                     inner_join(endpoint_mapping,by='ENDPOINT') %>%
#                     mutate(season_peak=case_when(month_peak_unadj<5 ~ 'winter',
#                                                  month_peak_unadj<9 ~ 'summer',
#                                                  TRUE ~ 'winter')) %>%
#                     filter(!is.na(brief)) %>%
#                     group_by(brief) %>%
#                     filter(n()>20) %>%
#                     ungroup() %>%
#                     count(brief,season_peak,name = 'count') %>%
#                     left_join(expand_grid(brief=unique(.$brief),season_peak=c('summer','winter')),.,by=c('brief','season_peak')) %>%
#                     mutate(count=ifelse(is.na(count),0,count)) %>%
#                     group_by(brief) %>%
#                     mutate(pct=100*count/sum(count)) %>%
#                     ungroup() %>%
#                     mutate(brief=ifelse(nchar(brief)>20,split_brief(brief),brief))
# filter(season_peak_dat,grepl('Infections|Pregnancy',brief))
#
#
# season_peak_plot <- ggplot(season_peak_dat,aes(brief,y=count,fill=season_peak)) +
#                     geom_col(position='dodge',alpha=0.5) +
#                     scale_y_continuous(expand=c(0,0)) +
#                     scale_fill_manual(values=c("#E18727FF","#0072B5FF"),name='Season peak') +
#                     xlab('') +
#                     ylab('Nr of endpoints') +
#                     theme_classic() +
#                     theme(axis.text.x = element_text(size=8,angle=-45,hjust=0),
#                           legend.position=c(0.8,0.8),
#                           plot.margin = margin(0,3,0,0, "cm"))


PTR_dat <- filter(seasonality_gam_summary_table,seasonality_pval_unadj<pval_thres) %>%
  pivot_longer(cols=c('ptr_est_unadj','ptr_est_adj'),names_to='ptr_type',values_to='estimate') %>%
  select(ENDPOINT,ptr_type,estimate) %>%
  inner_join(endpoint_mapping,by='ENDPOINT') %>%
  filter(!is.na(brief)) %>%
  group_by(brief) %>%
  filter(n()>40) %>%
  ungroup() %>%
  mutate(brief=ifelse(nchar(brief)>20,split_brief(brief),brief),
         ptr_type=factor(ptr_type,levels=c('ptr_est_unadj','ptr_est_adj'),labels=c('Unadjsted model','Adjusted model')))
  #filter(abs(delta_PTR)<1) %>%
  #reorder_by_median('delta_PTR',order='desc')

category_order <- filter(PTR_dat,ptr_type=='Adjusted model') %>%
                  group_by(brief) %>%
                  summarise(med=median(estimate)) %>%
                  ungroup() %>%
                  arrange(desc(med)) %>%
                  .$brief

PTR_plot <- ggplot(PTR_dat,aes(factor(brief,levels=category_order),y=estimate,fill=ptr_type)) +
            geom_boxplot(alpha=0.5) +
            geom_hline(yintercept=1,color='black',linetype='dashed') +
            coord_cartesian(ylim = c(1,4)) +
            scale_y_continuous(breaks=seq(1,4,by=0.5),expand=c(0,0.02)) +
            scale_fill_manual(values=c('#BC3C29FF','#0072B5FF'),name='') +
            xlab('') +
            ylab('Peak-to-trough ratio (PTR)') +
            theme_minimal() +
            theme(axis.text.x = element_text(angle=-45,hjust=0),
                  plot.margin = margin(0,0.5,0,0, "cm"),
                  legend.position=c(0.85,0.6))

#Delta PTR plots
# PTR_percentage_dat <- filter(seasonality_gam_summary_table,seasonality_pval_unadj<pval_thres) %>%
#   mutate(delta_PTR=-100*((ptr_est_unadj-ptr_est_adj)/(ptr_est_unadj-1))) %>%
#   select(ENDPOINT,delta_PTR) %>%
#   inner_join(endpoint_mapping,by='ENDPOINT') %>%
#   filter(!is.na(brief)) %>%
#   group_by(brief) %>%
#   filter(n()>20) %>%
#   ungroup() %>%
#   mutate(brief=ifelse(nchar(brief)>20,split_brief(brief),brief)) %>%
#   #filter(abs(delta_PTR)<1) %>%
#   reorder_by_median('delta_PTR',order='desc')
#
# PTR_percentage_plot <- group_by(PTR_percentage_dat,brief) %>%
#   summarise(median(delta_PTR)) %>%
#   #slice(1:7) %>%
#   inner_join(PTR_percentage_dat,by='brief') %>%
#   ggplot(aes(brief,y=delta_PTR)) +
#   geom_boxplot(alpha=0.5,fill="#0072B5FF") +
#   geom_hline(yintercept=0,color='red') +
#   #scale_y_continuous(limits=c(0,NA),expand=c(0.01,0)) +
#   xlab('') +
#   ylab('Delta PTR \\%') +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=-45,hjust=0),
#         plot.margin = margin(0,3,0,0, "cm"))

heartfailure_plot <- autoplot(seasonality_gam_FinRegistry_unadj[['I9_HEARTFAIL_ALLCAUSE']],type='seasonality') + ggtitle('')

set.seed(10)
dummy_boxplot <- tibble(Norway=rnorm(10),Denmark=rnorm(10),Iceland=rnorm(10)) %>%
                 pivot_longer(cols = everything()) %>%
                 ggplot(aes(x=name,y=value,col=name)) +
                 geom_boxplot(width=0.6,linewidth=1) +
                 scale_color_manual(values=c('#BC3C29FF','#0072B5FF','#E18727FF')) +
                 xlab('') +
                 ylab('Timing difference of peak') +
                 theme_minimal() +
                 theme(#axis.text.x = element_text(size=12),
                        legend.position = 'none')

seasonal_overview_plot <- ((avg_spline_plot  + ggtitle('A'))  + (peak_trough_hist + ggtitle('B'))) /
                          ((PTR_plot + ggtitle('C'))) /
                          ((dummy_boxplot + ggtitle('D'))+(heartfailure_plot + ggtitle('E')))
ggsave(filename='seasonality_paper/results/png/seasonality_overview.png',plot = seasonal_overview_plot,device='png',width = 7.2,height=7.5)
save_tikz_plot(seasonal_overview_plot,width=7.2,height=7.5,filename = 'seasonality_paper/results/tex/figure2_seasonal_overview_plot.tex')

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


#month peak
# month_peak_dat <- filter(seasonal_summary,(smooth_pval)<pval_thres) %>%
#   mutate(month_peak=ifelse(month_peak<1,floor(month_peak+12),floor(month_peak))) %>%
#   count(month_peak) %>%
#   mutate(pct=100*n/sum(n))
#
# filter(month_peak_dat,month_peak>=4 & month_peak<=8) %>%
#   summarise(sum(pct))
#
# filter(month_peak_dat,month_peak %in% c(2,3,10)) %>%
#   summarise(sum(pct))
#
# month_peak_plot <- ggplot(month_peak_dat,aes(month_peak,y=pct)) +
#   geom_col(position='dodge',fill="#BC3C29FF",alpha=0.5,width=0.8) +
#   scale_x_continuous(breaks=seq(1,12),expand=c(0,0.2)) +
#   scale_y_continuous(expand=c(0,0)) +
#   xlab('Peak month number') +
#   ylab('\\% of endpoints') +
#   theme_classic() +
#   theme(plot.margin = margin(0,0,0,0, "cm"),
#         axis.title.x = element_text(vjust = 1))


