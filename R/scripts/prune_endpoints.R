library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)
source('seasonality_dashboard/utils/seasonality_utils.R')
source('R/scripts/utils.R')
monthly_counts_FinRegistry  <- read_tsv('data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
                               filter_monthly_counts(.,endpoint_mapping_path='data/finngen_R10_endpoint_core_noncore_1.0.xlsx')
seasonal_summary_FinRegistry <- read_tsv('data/FINREGISTRY_seasonal_summary.txt')
seasonal_summary_FinRegistry_adj <- read_tsv('data/FINREGISTRY_seasonal_summary_mean_offset_adj.txt')
cases_per_season_FinnGen <- read_tsv('data/FINNGEN_R11_cases_per_season.txt')
GWAS_summary_FinnGen <- read_tsv('data/finngen_R10_pheno_n.tsv')
ldsc_heritability_FinnGen <- read_tsv('data/finngen_R11_FIN.ldsc.heritability.tsv')


endpoint_categories <- read_excel('data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=2)
endpoint_table <- read_excel('data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=1,range = cell_cols("A:F")) %>%
                  separate_rows(TAGS,sep=',') %>%
                  mutate(TAGS=gsub(' +','',TAGS)) %>%
                  filter(TAGS!='#AVO') %>% # Remove avohilmo tag,as it is irrelevant to our study
                  left_join(endpoint_categories,by=c('TAGS'='code')) %>%
                  select(NAME,TAGS,LEVEL,LONGNAME,brief,CORE_ENDPOINTS)

endpoint_mapping  <- distinct(monthly_counts_FinRegistry,ENDPOINT) %>%
                     left_join(distinct(endpoint_table,ENDPOINT=NAME,brief,TAGS,LEVEL,CORE_ENDPOINTS,LONGNAME),by='ENDPOINT')


endpoint_brief_remove <-  filter(endpoint_categories,!is.na(brief)) %>%
                        filter(grepl('Katri',brief)) %>% .$brief
endpoint_brief_remove_strict <- c('External causes','Drug purchase','Drug-induced adverse effects','Operation',
                                'Health status','Pregnancy and childbirth','Musculoskeletal and connective tissue',
                                'Dental','Perinatal','Genitourinary system','Miscellaneous, not yet classified endpoints',
                                'Other, not yet classified endpoints (same as #MISC)','Common')


#filter(CORE_ENDPOINTS=='yes') %>%
filtered_endpoints <- inner_join(seasonal_summary_FinRegistry,select(seasonal_summary_FinRegistry_adj,ENDPOINT,smooth_pval_adj=smooth_pval),by='ENDPOINT') %>%
                      filter(smooth_pval<(0.05/n()),
                             smooth_pval_adj<(0.05/n())) %>%
                      inner_join(cases_per_season_FinnGen,by='ENDPOINT') %>%
                      filter(cases>10000,pmin(cases_high,cases_low)>4000) %>%
                      inner_join(endpoint_mapping,by='ENDPOINT') %>%
                      filter(!(brief %in% endpoint_brief_remove)) %>%
                      group_by(ENDPOINT) %>%
                      filter(!any(brief %in% endpoint_brief_remove_strict)) %>%
                      filter(if(n()>1) c(T,rep(F,n()-1)) else T) %>%
                      ungroup() %>%
                      inner_join(select(GWAS_summary_FinnGen,ENDPOINT=phenocode,num_gw_significant),by='ENDPOINT') %>%
                      inner_join(select(ldsc_heritability_FinnGen,ENDPOINT=PHENO,H2),by='ENDPOINT') %>%
                      filter(num_gw_significant>0)


#Genetic correlation pruning
genetic_cor_dat <- read_tsv('data/finngen_R11_FIN.ldsc.summary.tsv') %>%
                   filter(p1 %in% filtered_endpoints$ENDPOINT & p2 %in% filtered_endpoints$ENDPOINT)
endpoints_left <- filtered_endpoints$ENDPOINT
endpoint_clusters <- list()
i=1
genetic_cor_dat_endpoints <- genetic_cor_dat
while(i<=length(endpoints_left)){
  correlated_endpoints <- filter(genetic_cor_dat_endpoints,p1==endpoints_left[i] | p2==endpoints_left[i]) %>%
                          filter(CONVERGED) %>%
                          filter(p<0.05,abs(rg)>0.9) %>%
                          pivot_longer(cols=c('p1','p2'),values_to='ENDPOINT') %>%
                          filter(ENDPOINT!=endpoints_left[i])
  endpoints_left <- endpoints_left[!(endpoints_left %in% correlated_endpoints$ENDPOINT)]
  genetic_cor_dat_endpoints <- filter(genetic_cor_dat,p1 %in% endpoints_left & p2 %in% endpoints_left)
  endpoint_clusters[[endpoints_left[i]]] <- correlated_endpoints$ENDPOINT
  i=i+1
}
selected_endpoints <- tibble(ENDPOINT=endpoints_left) %>%
                      inner_join(endpoint_mapping,by='ENDPOINT') %>%
                      group_by(ENDPOINT) %>%
                      summarise(ENDPOINT_category=paste(brief,collapse=', ')) %>%
                      inner_join(filtered_endpoints,.,by='ENDPOINT') %>%
                      mutate(cluster_endpoints=sapply(endpoint_clusters,paste,collapse=', ')) %>%
                      select(ENDPOINT,cases,cases_low,cases_high,ptr_est,smooth_pval,smooth_pval_adj,H2,num_gw_significant,cluster_endpoints,LONGNAME,ENDPOINT_category)

write.table(selected_endpoints,file='results/seasonality_GWAS_endpoints.tsv',quote=F,row.names=F,sep='\t')
write_xlsx(selected_endpoints,path='results/seasonality_GWAS_endpoints.xlsx')
