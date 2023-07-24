library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
source('seasonality_dashboard/utils/seasonality_utils.R')
source('R/scripts/utils.R')
freeze <- 'R11'
pheno_dat <- read_tsv(paste0('~/Documents/thesis/data/FINNGEN_',freeze,'_endpoints_individual_dates.txt'))
seasonal_smooth_terms_FinRegistry <- read_tsv('~/Documents/thesis/data/FINREGISTRY_seasonal_splines.txt')

cases_per_season <- lapply(unique(pheno_dat$ENDPOINT),function(e){
  if(e %in% seasonal_smooth_terms_FinRegistry$ENDPOINT){
    seasonal_pheno <- get_seasonal_phenotype(endpoint_id=e,pheno_dat=pheno_dat,seasonal_spline=seasonal_smooth_terms_FinRegistry,spline_endpoint = e)
    tibble(ENDPOINT=e,cases=nrow(seasonal_pheno),cases_high=sum(seasonal_pheno$seasonal_val_binary==1),cases_low=sum(seasonal_pheno$seasonal_val_binary==0))
  }else{
    return(NULL)
  }
}) %>% bind_rows()

write.table(cases_per_season,paste0('data/FINNGEN_',freeze,'_cases_per_season.txt'),row.names=F,sep='\t',quote=F)
