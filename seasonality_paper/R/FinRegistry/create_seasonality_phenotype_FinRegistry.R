.libPaths('/shared-directory/sd-tools/apps/R/lib/')
library(dplyr)
library(readr)
library(data.table)
library(gseasonality)
library(lubridate)

batch_read_file = function(path,batch_size=5e5,cols){
  con_path=file(path,'rb')
  col_names = gsub('\"','',unlist(strsplit(readLines(con_path,n=1),split=',')))
  i=1
  dat=list()
  while(TRUE){
    print(i)
    dat_batch = fread(text=readLines(con_path,n = batch_size),sep=',')
    if(nrow(dat_batch)==0) break
    names(dat_batch) = col_names
    dat[[i]] = dat_batch[,..cols]
    i=i+1
  }
  close(con_path)
  return(rbindlist(dat))
}

year_start = 1999
year_end = 2019

seasonality_gam_FinRegistry_unadj = readRDS(paste0('results/seasonality_gam_FinRegistry_unadj_',year_start,'-',year_end,'.rds'))
pval_thres=0.05/length(seasonality_gam_FinRegistry_unadj)
for(e in names(seasonality_gam_FinRegistry_unadj)){
  if(seasonality_gam_FinRegistry_unadj[[e]]$summary$pval>=pval_thres){
    seasonality_gam_FinRegistry_unadj[[e]] = NULL
  }
}
starting_letter_distr=table(substr(names(seasonality_gam_FinRegistry_unadj),1,1))
starting_letters=sort(names(starting_letter_distr))

small_groups=starting_letter_distr<20
starting_letter_regex = c(paste0('^',starting_letters[!small_groups]),paste(paste0('^',starting_letters[small_groups],collapse='|')))

phenotype_minimal = batch_read_file(path='~/Projects/SD-Connect/project_2007099/minimal_phenotype_20230110/minimal_phenotype_20221216.csv',
                                    batch_size=5e5,cols=c('FINREGISTRYID','date_of_birth','sex','death_date'))


for(i in 1:length(starting_letter_regex)){
  print(starting_letter_regex[i])
  system(paste0('awk \'BEGIN{FS=","}{if(NR==1 || $2~/',starting_letter_regex[i],
                '/){print}}\' ~/Projects/SD-Connect/project_2007099/endpointer_20221112/densified_first_events_DF10_no_omits_2022-09-20.txt > ',tempdir(),'/first_event_tmp.txt'))
  first_endpoint_dat = batch_read_file(paste0(tempdir(),'/first_event_tmp.txt'),batch_size = 1e6,cols=c('FINREGISTRYID','ENDPOINT','AGE'))
  endpoints = sort(unique(first_endpoint_dat$ENDPOINT))
  endpoints = endpoints[endpoints %in% names(seasonality_gam_FinRegistry_unadj)]
  seasonality_phenotype_FinRegistry = lapply(endpoints,function(endpoint){
    print(endpoint)
    endpoint_dates = first_endpoint_dat[ENDPOINT==endpoint,]
    dates_dat = phenotype_minimal[endpoint_dates,on='FINREGISTRYID',nomatch=0]
    idx_leap_date = which(grepl('^29-02-',dates_dat$date_of_birth))
    idx_invalid_leap_date = idx_leap_date[which(is.na(dmy(dates_dat$date_of_birth[idx_leap_date])))]
    substring(dates_dat$date_of_birth[idx_invalid_leap_date],1,2)='28' #invalid dates with 29-02 in non leap years: set to 28-02
    dates_dat[,EVENT_DATE:=as_date(dmy(date_of_birth) + dyears(AGE))]
    dates_dat = dates_dat[,.(ID=FINREGISTRYID,ENDPOINT,EVENT_DATE,AGE,sex,date_of_birth,death_date)]
    seasonality_pheno = as.data.table(get_seasonality_phenotype(mod=seasonality_gam_FinRegistry_unadj[[endpoint]],
                                                  data=as_tibble(dates_dat),year_start=1999,year_end=2019))
    dates_dat = dates_dat[seasonality_pheno,on=c('ID','EVENT_DATE')]
    dates_dat[,n:=nrow(dates_dat)]
  }) %>% rbindlist()
  fwrite(seasonality_phenotype_FinRegistry[,.(ID,ENDPOINT,EVENT_DATE,AGE,sex,date_of_birth,death_date,pheno_qt,pheno_binary)],
         file = paste0('results/seasonality_phenotype_FinRegistry_',year_start,'-',year_end,'.tsv'),sep="\t",quote=F,append = i!=1,col.names = i==1)
  fwrite(unique(seasonality_phenotype_FinRegistry,by=c('ENDPOINT','n'))[,.(ENDPOINT,n)],
         file=paste0('results/seasonality_phenotype_size_FinRegistry_',year_start,'-',year_end,'.tsv'),sep='\t',quote=F,append=i!=1,col.names = i==1)
}
