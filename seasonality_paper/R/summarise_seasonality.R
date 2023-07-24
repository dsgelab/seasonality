library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(gseasonality)
source('seasonality_paper/R/utils.R')
year_start <- 1999
year_end <- 2019
monthly_counts_FinRegistry  <- read_tsv('seasonality_paper/data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
  filter_monthly_counts(.,year_start=year_start,year_end=year_end)

endpoints <- unique(monthly_counts_FinRegistry$ENDPOINT)
mod_list_unadj <- lapply(endpoints,function(e){
  print(e)
  seasonality_gam(monthly_counts=filter(monthly_counts_FinRegistry,ENDPOINT==e),year_start=year_start,year_end=year_end,adjusted = F)
})
names(mod_list_unadj) <- endpoints
prepare_output(mod_list_unadj,path=paste0('seasonality_paper/results/seasonality_gam_FinRegistry_unadj_',year_start,'-',year_end,'.rds'))

mod_list_adj <- lapply(endpoints,function(e){
  print(e)
  seasonality_gam(monthly_counts=filter(monthly_counts_FinRegistry,ENDPOINT==e),year_start=year_start,year_end=year_end,adjusted = T)
})
names(mod_list_adj) <- endpoints
prepare_output(mod_list_adj,path=paste0('seasonality_paper/results/seasonality_gam_FinRegistry_adj_',year_start,'-',year_end,'.rds'))
