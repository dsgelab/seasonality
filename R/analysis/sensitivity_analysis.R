library(readr)
source('seasonality_dashboard/utils/seasonality_utils.R')
GWAS_info_dat <- read_tsv('data/finngen_R10_finngen_R10_analysis_data_finngen_R10_pheno_n.tsv')

monthly_counts_FinRegistry  <- read_tsv('data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
                                 filter_monthly_counts(.,GWAS_info_dat=GWAS_info_dat,endpoint_mapping_path='data/finngen_R10_endpoint_core_noncore_1.0.xlsx')


a <- 1;b <- 13
###### Seasonal splines for all sensitivity analyses
seasonal_splines_after2000 <- filter(monthly_counts_FinRegistry,EVENT_YEAR>=2000) %>%
                    group_by(ENDPOINT) %>%
                    group_modify(~extract_seasonal_spline(.x,a=a,b=b))
write.table(seasonal_splines_after2000,file='data/FINREGISTRY_seasonal_splines_after2000.txt',sep='\t',row.names=F,quote=F)



seasonal_splines_k4 <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                          group_modify(~extract_seasonal_spline(.x,a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=4))
write.table(seasonal_splines_k4,file='data/FINREGISTRY_seasonal_splines_k4.txt',sep='\t',row.names=F,quote=F)



seasonal_splines_k8 <- group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
  group_modify(~extract_seasonal_spline(.x,a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=8))
write.table(seasonal_splines_k8,file='data/FINREGISTRY_seasonal_splines_k8.txt',sep='\t',row.names=F,quote=F)

###### Seasonality summary for all sensitivity analyses

seasonal_summary_FinRegistry_after2000 = filter(monthly_counts_FinRegistry,EVENT_YEAR>=2000) %>%
                                         group_by(ENDPOINT) %>%
                                         group_modify(~summarise_seasonality(.x,a=a,b=b,adjustment='',seasonal_spline_avg=NULL)) %>%
                                         arrange(log10_ptr_pval)

write.table(seasonal_summary_FinRegistry_after2000,file='data/FINREGISTRY_seasonal_summary_after2000.txt',sep='\t',row.names=F,quote=F)


seasonal_summary_FinRegistry_k4 = group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                                  group_modify(~summarise_seasonality(.x,a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=4,seasonal_spline_avg=NULL)) %>%
                                  arrange(log10_ptr_pval)

write.table(seasonal_summary_FinRegistry_k4,file='data/FINREGISTRY_seasonal_summary_k4.txt',sep='\t',row.names=F,quote=F)


seasonal_summary_FinRegistry_k8 = group_by(monthly_counts_FinRegistry,ENDPOINT) %>%
                                         group_modify(~summarise_seasonality(.x,a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=8,seasonal_spline_avg=NULL)) %>%
                                         arrange(log10_ptr_pval)

write.table(seasonal_summary_FinRegistry_k8,file='data/FINREGISTRY_seasonal_summary_k8.txt',sep='\t',row.names=F,quote=F)
