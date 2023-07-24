library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)
library(gseasonality)
source('seasonality_paper/R/utils.R')

freeze <- 'R11'
year_start <- 1999
year_end <- 2019
monthly_counts_FinRegistry  <- read_tsv('seasonality_paper/data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
                               filter_monthly_counts(.,year_start=year_start,year_end=year_end)
seasonality_FinRegistry_unadj <- readRDS(paste0('seasonality_paper/results/seasonality_gam_FinRegistry_unadj_',year_start,'-',year_end,'.rds'))
seasonality_FinRegistry_adj <- readRDS(paste0('seasonality_paper/results/seasonality_gam_FinRegistry_adj_',year_start,'-',year_end,'.rds'))

# endpoint_dates_FinnGen <- read_tsv(paste0('seasonality_paper/data/FINNGEN_',freeze,'_endpoints_individual_dates.txt'))
seasonality_phenotypes_FinnGen <- read_tsv(paste0('seasonality_paper/results/seasonality_phenotype_FinnGen_',freeze,'_unadj_',year_start,'-',year_end,'.tsv'))
GWAS_summary_FinnGen <- read_tsv(paste0('seasonality_paper/data/finngen_',freeze,'_pheno_n.tsv'))
ldsc_heritability_FinnGen <- read_tsv(paste0('seasonality_paper/data/finngen_',freeze,'_FIN.ldsc.heritability.tsv'))

endpoint_categories <- read_excel('seasonality_paper/data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=2)
endpoint_mapping_manual <- read_excel('seasonality_paper/data/finngen_R10_endpoint_core_noncore_1.0.xlsx',sheet=1) %>%
                          select(NAME,TAGS,LONGNAME) %>%
                          separate_rows(TAGS,sep=',') %>%
                          left_join(endpoint_categories,by=c('TAGS'='code')) %>%
                          filter(TAGS!='#AVO') %>% # Remove avohilmo tag,as it is irrelevant to our study
                          select(NAME,TAGS,LONGNAME,brief)
endpoint_mapping_prefix <- get_endpoint_mapping(names(seasonality_gam_FinRegistry_unadj),endpoint_categories)
endpoint_mapping  <- left_join(endpoint_mapping_prefix,endpoint_mapping_manual,by=c('ENDPOINT'='NAME'),suffix=c('_prefix','_manual')) %>%
                     mutate(brief_prefix=ifelse(is.na(brief_prefix),brief_manual,brief_prefix),
                            brief_manual=ifelse(is.na(brief_manual),LONGNAME,brief_manual),
                            brief=ifelse(brief_prefix!=brief_manual,brief_prefix,brief_manual)) %>%
                     distinct(ENDPOINT,brief,prefix)



endpoint_brief_remove <-  filter(endpoint_categories,!is.na(brief)) %>%
                          filter(grepl('Katri',brief)) %>% .$brief
endpoint_brief_remove_strict <- c('External causes','Drug purchase','Drug-induced adverse effects','Operation',
                                  'Health status','Pregnancy and childbirth','Musculoskeletal and connective tissue','Endocrine, nutritional and metabolic',
                                  'Dental','Perinatal','Genitourinary system','Miscellaneous, not yet classified endpoints',
                                  'Other, not yet classified endpoints (same as #MISC)','Common')

seasonality_summary <- lapply(names(seasonality_FinRegistry_unadj),function(e){
                          tibble(ENDPOINT=e,
                                 ptr_est=seasonality_FinRegistry_unadj[[e]]$summary$ptr$estimate,
                                 seasonality_pval=seasonality_FinRegistry_unadj[[e]]$summary$pval,
                                 seasonality_pval_adj=seasonality_FinRegistry_adj[[e]]$summary$pval)
                        }) %>% bind_rows()

cases_per_season_FinnGen <- count(seasonality_phenotypes_FinnGen,ENDPOINT,pheno_binary) %>%
                            pivot_wider(id_cols='ENDPOINT',names_from='pheno_binary',values_from='n',names_prefix = 'cases_')

filtered_endpoints <- filter(seasonality_summary,seasonality_pval<(0.05/n()),
                                                 seasonality_pval_adj<(0.05/n())) %>%
                      inner_join(cases_per_season_FinnGen,by='ENDPOINT') %>%
                      filter(cases>10000,pmin(cases_1,cases_0)>4000) %>%
                      inner_join(endpoint_mapping,by='ENDPOINT') %>%
                      filter(!(brief %in% endpoint_brief_remove)) %>%
                      filter(!grepl('REIMB',ENDPOINT)) %>% #remove reimbursement endpoints not caught by previous filtering
                      group_by(ENDPOINT) %>%
                      filter(!any(brief %in% endpoint_brief_remove_strict)) %>%
                      filter(if(n()>1) c(T,rep(F,n()-1)) else T) %>%
                      ungroup() %>%
                      inner_join(select(GWAS_summary_FinnGen,ENDPOINT=phenocode,num_gw_significant),by='ENDPOINT') %>%
                      inner_join(select(ldsc_heritability_FinnGen,ENDPOINT=PHENO,H2),by='ENDPOINT') %>%
                      filter(num_gw_significant>0) %>%
                      arrange(desc(ptr_est))


#Genetic correlation pruning
genetic_cor_dat <- read_tsv('seasonality_paper/data/finngen_R11_FIN.ldsc.summary.tsv') %>%
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
  select(ENDPOINT,cases,cases_low,cases_high,ptr_est,seasonality_pval,seasonality_pval_adj,H2,num_gw_significant,cluster_endpoints,LONGNAME,ENDPOINT_category)

write.table(selected_endpoints,file='results/seasonality_GWAS_endpoints.tsv',quote=F,row.names=F,sep='\t')
write_xlsx(selected_endpoints,path='results/seasonality_GWAS_endpoints.xlsx')
