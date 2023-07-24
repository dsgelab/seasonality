library(readr)
source('seasonality_dashboard/utils/seasonality_utils.R')
monthly_counts_FinRegistry  <- read_tsv('data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
                               filter_monthly_counts(.,endpoint_mapping_path='data/finngen_R10_endpoint_core_noncore_1.0.xlsx')

a <- 1;b <- 13

####### Extract seasonal splines ########
seasonal_splines <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                    group_modify(~extract_seasonal_spline(.x,a=a,b=b))
write.table(seasonal_splines,file='data/FINREGISTRY_seasonal_splines.txt',sep='\t',row.names=F,quote=F)
write.table(seasonal_splines,file='seasonality_dashboard/data/FINREGISTRY_seasonal_splines.txt',sep='\t',row.names=F,quote=F)

seasonal_splines_binary_adj <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
  group_modify(~extract_seasonal_spline(.x,a=a,b=b,adjustment='binary'))
write.table(seasonal_splines_binary_adj,file='data/FINREGISTRY_seasonal_splines_binary_adj.txt',sep='\t',row.names=F,quote=F)

#Month length adjusted
seasonal_splines_month_adj <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                    group_modify(~extract_seasonal_spline(.x,a=a,b=b,adj_month_length = T))
write.table(seasonal_splines_month_adj,file='data/FINREGISTRY_seasonal_splines_month_adj.txt',sep='\t',row.names=F,quote=F)


# Extract the average seasonal smooth
seasonal_spline_avg <- group_by(seasonal_splines,month) %>%
  summarise(avg_seasonal_val=mean(seasonal_val))

######### Seasonal summary ############
## Unadjusted
seasonal_summary <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                    group_modify(~summarise_seasonality(.x,a=a,b=b,mod_null='noadj')) %>%
                    arrange(desc(ptr_est))
write.table(seasonal_summary,file='data/FINREGISTRY_seasonal_summary.txt',sep='\t',row.names=F,quote=F)
write.table(seasonal_summary,file='seasonality_dashboard/data/FINREGISTRY_seasonal_summary.txt',sep='\t',row.names=F,quote=F)

## Month adjusted
seasonal_summary_month_adj <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                              group_modify(~summarise_seasonality(.x,a=a,b=b,adj_month_length = T,mod_null='noadj')) %>%
                              arrange(desc(ptr_est))
write.table(seasonal_summary_month_adj,file='data/FINREGISTRY_seasonal_summary.txt',sep='\t',row.names=F,quote=F)

## Binary adjusted
seasonal_summary_binary_adj = group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                              group_modify(~summarise_seasonality(.x,a=a,b=b,adjustment='binary',mod_null='noadj',seasonal_spline_avg=seasonal_spline_avg)) %>%
                              arrange(desc(ptr_est))

write.table(seasonal_summary_binary_adj,file='data/FINREGISTRY_seasonal_summary_binary_adj.txt',sep='\t',row.names=F,quote=F)

## Mean offset adjusted
seasonal_summary_mean_offset_adj = group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                                   group_modify(~summarise_seasonality(.x,a=a,b=b,adjustment='mean_offset',mod_null='noadj',seasonal_spline_avg=seasonal_spline_avg)) %>%
                                   arrange(desc(ptr_est))

write.table(seasonal_summary_mean_offset_adj,file='data/FINREGISTRY_seasonal_summary_mean_offset_adj.txt',sep='\t',row.names=F,quote=F)

## Mean covariate adjusted
seasonal_summary_mean_covariate_adj = group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                                      group_modify(~summarise_seasonality(.x,a=a,b=b,adjustment='mean_covariate',mod_null='noadj',seasonal_spline_avg=seasonal_spline_avg)) %>%
                                      arrange(desc(ptr_est))

write.table(seasonal_summary_mean_covariate_adj,file='data/FINREGISTRY_seasonal_summary_mean_covariate_adj.txt',sep='\t',row.names=F,quote=F)



