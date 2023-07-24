library(dplyr)
library(readr)
library(gseasonality)
freeze <- 'R11'
year_start <- 1999
year_end <- 2019
endpoint_dates_FinnGen <- read_tsv(paste0('seasonality_paper/data/FINNGEN_',freeze,'_endpoints_individual_dates.txt'))

seasonality_gam_FinRegistry_unadj <- readRDS(paste0('seasonality_paper/results/seasonality_gam_FinRegistry_unadj_',year_start,'-',year_end,'.rds'))
seasonality_gam_FinRegistry_adj <- readRDS(paste0('seasonality_paper/results/seasonality_gam_FinRegistry_adj_',year_start,'-',year_end,'.rds'))

endpoints_FinRegistry <- names(seasonality_gam_FinRegistry_unadj)
endpoints_FinnGen <- unique(endpoint_dates_FinnGen$ENDPOINT)
seasonality_phenotype_unadj_FinnGen <- lapply(endpoints_FinRegistry,function(e){
  if(e %in% endpoints_FinnGen){
    seasonality_pheno <- get_seasonality_phenotype(mod=seasonality_gam_FinRegistry_unadj[[e]],
                                                   data=filter(endpoint_dates_FinnGen,ENDPOINT==e) %>%
                                                        select(ID=FINNGENID,EVENT_DATE)) %>%
                         mutate(ENDPOINT=e) %>%
                         select(ENDPOINT,everything())
    return(seasonality_pheno)
  }else{
    return(NULL)
  }
}) %>% bind_rows()

seasonality_phenotype_adj_FinnGen <- lapply(endpoints_FinRegistry,function(e){
  if(e %in% endpoints_FinnGen){
    seasonality_pheno <- get_seasonality_phenotype(mod=seasonality_gam_FinRegistry_adj[[e]],
                                                   data=filter(endpoint_dates_FinnGen,ENDPOINT==e) %>%
                                                     select(ID=FINNGENID,EVENT_DATE)) %>%
                         mutate(ENDPOINT=e) %>%
                         select(ENDPOINT,everything())
    return(seasonality_pheno)
  }else{
    return(NULL)
  }
}) %>% bind_rows()

seasonality_phenotype_FinnGen <- inner_join(seasonality_phenotype_unadj_FinnGen,seasonality_phenotype_adj_FinnGen,by=c('ID','EVENT_DATE','EVENT_MONTH_DEC'),suffix=c('','_adj'))
rm(seasonality_phenotype_adj_FinnGen,
   seasonality_phenotype_unadj_FinnGen,
   seasonality_gam_FinRegistry_unadj,
   seasonality_gam_FinRegistry_adj)
write_tsv(seasonality_phenotype_FinnGen,
          file = paste0('seasonality_paper/results/seasonality_phenotype_FinnGen_',freeze,'_',year_start,'-',year_end,'.tsv'))
