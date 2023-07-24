.libPaths('/shared-directory/sd-tools/apps/R/lib/')
library(data.table)
library(comorbidity)
library(dplyr)

seasonality_gam_FinRegistry = readRDS('results/seasonality_gam_FinRegistry_unadj_1999-2019.rds')
endpoints = lapply(names(seasonality_gam_FinRegistry), function(n){
  if(seasonality_gam_FinRegistry[[n]]$summary$pval<0.05/length(seasonality_gam_FinRegistry)){
    return(data.table(ENDPOINT=n))
  }else{
    return(NULL)
  }
}) %>% rbindlist()

comorbs=c('ami','chf','pvd','cevd','dementia','copd','rheumd','pud','mld',
          'diab','diabwc','hp','rend','canc','msld','metacanc','aids')
weights=as.matrix(c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,6,6))
ICD_versions = c('9','10')

detailed_path='Projects/SD-Connect/project_2007099/detailed_longitudinal_20221118/detailed_longitudinal_DF10_2022-11-11.csv'
con_detailed=file(detailed_path,'rb')
col_names_detailed = unlist(strsplit(readLines(con_detailed,n=1),split=','))

CCI_dat = setDT(NULL)
next_iter_dat = data.table(FINREGISTRYID=character(),EVENT_AGE=numeric(),PRIMARY_ICD=character(),ICD_VERSION=character())
batch_size = 1e6
curr_batch = 1
t = Sys.time()
while(TRUE){
  print(curr_batch)
  dat_detailed = fread(text=readLines(con_detailed,n = batch_size),sep=',')
  if(nrow(dat_detailed)==0) break
  names(dat_detailed) = col_names_detailed
  dat_detailed[,ICDver:=ifelse(ICDVER=='03','10',ICDVER)]
  dat_detailed = dat_detailed[SOURCE!='PURCH' & ICDVER %in% ICD_versions,.(FINREGISTRYID,EVENT_AGE,PRIMARY_ICD=CODE1,ICD_VERSION=ICDVER)]
  dat_detailed = rbind(next_iter_dat,dat_detailed)
  next_iter_dat = dat_detailed[FINREGISTRYID==tail(FINREGISTRYID,1),]
  dat_detailed = dat_detailed[!next_iter_dat,on=c('FINREGISTRYID','EVENT_AGE','PRIMARY_ICD','ICD_VERSION')]
  setorder(dat_detailed,FINREGISTRYID,EVENT_AGE)
  dat_detailed = unique(dat_detailed,by = c('FINREGISTRYID','PRIMARY_ICD'))
  dat_detailed[,pseudo_ID:=paste(FINREGISTRYID,EVENT_AGE,sep='_')]
  comorb_dat = lapply(sort(unique(dat_detailed$ICD_VERSION)),function(version){
    as.data.table(comorbidity::comorbidity(dat_detailed[ICD_VERSION==version,],id='pseudo_ID',code='PRIMARY_ICD',map=paste0('charlson_icd',version,'_quan'),assign0=F))
  }) %>% rbindlist()
  comorb_dat[,FINREGISTRYID:=gsub('_.*','',pseudo_ID)]
  comorb_dat[,EVENT_AGE:=as.numeric(gsub('.*_','',pseudo_ID))]
  setorder(comorb_dat,FINREGISTRYID,EVENT_AGE)
  comorb_dat_cum = comorb_dat[,lapply(.SD,cummax),by=FINREGISTRYID,.SD=comorbs]
  CCI_dat_batch = data.table(FINREGISTRYID=comorb_dat_cum$FINREGISTRYID,
                       AGE=as.numeric(sub('.*_','',comorb_dat$pseudo_ID)),
                       CCI=as.numeric(as.matrix(comorb_dat_cum[,..comorbs]) %*% weights))
  CCI_dat_batch = unique(CCI_dat_batch,by=c('FINREGISTRYID','CCI'))
  CCI_dat = rbind(CCI_dat,CCI_dat_batch)
  curr_batch = curr_batch + 1
}
Sys.time()-t
close(con_detailed)
fwrite(CCI_dat,file = 'results/CCI_progression_FinRegistry.tsv',sep="\t",quote=F)


