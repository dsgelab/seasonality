#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(lubridate)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("-e", "--endpoint"), type="character", default='I9_HEARTFAIL_ALLCAUSE',
              help="Which endpoint should be analysed", metavar="character"),
  make_option(c("-t", "--type"), type="character", default="co_binary",
              help=paste0("Which type of GWAS to perform. Available types are case-only",
                          "quantitative trait ('co_qt'), case-only binary trait ('co_binary'), ",
                          "case-control binary trait for highs (cc_binary_high), ",
                          "case-control binary trait for lows (cc_binary_low),
                            and a normal case-control (cc_binary_all)"), metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

get_seasonal_phenotype <- function(endpoint_id,pheno_dat,seasonal_spline,spline_endpoint=endpoint_id){
  dates <- seq(ymd('1998-01-01'),ymd('2019-12-31'),by='1 day')
  dates_dat <- tibble(year=year(dates),month=month(dates),day=day(dates)) %>%
               group_by(year,month) %>%
               summarise(num_days=max(day))
  pheno_dat_endpoint <- as_tibble(pheno_dat[ENDPOINT==endpoint_id,]) %>%
                        rename(FID=FINNGENID) %>%
                        mutate(EVENT_MONTH=month(EVENT_DATE),
                               EVENT_DAY=day(EVENT_DATE)) %>%
                        inner_join(dates_dat,by=c('EVENT_YEAR'='year','EVENT_MONTH'='month')) %>%
                        mutate(EVENT_MONTH_DEC=EVENT_MONTH + EVENT_DAY/num_days) %>%
                        mutate(EVENT_MONTH_DEC=ifelse(EVENT_MONTH_DEC>12.5,
                                                      EVENT_MONTH_DEC-12,
                                                      EVENT_MONTH_DEC))
  seasonal_spline_endpoint <- filter(seasonal_spline,ENDPOINT==spline_endpoint)
  pheno_dat_endpoint <- mutate(pheno_dat_endpoint,
                               seasonal_val=approx(x=seasonal_spline_endpoint$month,
                                                   y=seasonal_spline_endpoint$seasonal_val,
                                                   xout=EVENT_MONTH_DEC)$y,
                               seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                               seasonal_val_binary=as.integer(seasonal_val > 0))
  return(pheno_dat_endpoint)
}
#Create folder structure for GWAS
base_path <- '/finngen/red/'
user_path <- 'solviro/GWAS_R10/REGENIE/'
analysis_id <- paste(args$endpoint,args$type,sep='.')
GWAS_path_tmp <- paste0('~/Documents/thesis/tmp/')
template_path <- paste0(base_path,user_path,'template/')
phenolist_filename <-  paste0('pheno_',analysis_id,'.txt')
pheno_filename <-  paste0(analysis_id,'_R10_GWAS.txt.gz')
#create directory of not yet created
if(!dir.exists(GWAS_path)){
  dir.create(GWAS_path)
}
#Add infrastructure for REGENIE
#copy regenie files from template directory
files <- c('finngen_r10_bgen_sample_list.txt','r10_bgen_index_list.txt',
           'regenie.wdl','regenie_sub.zip')
for(f in files){
  file.copy(from=paste0(template_path,f),to=GWAS_path_tmp,overwrite = TRUE)
}
#Edit json file inline with the current GWAS to be run
json_lines <- readLines(paste0(template_path,'regenie_ENDPOINT.json')) %>%
              gsub('USER_PATH',paste0(user_path,analysis_id,'/'),.) %>%
              gsub('PHENOLIST_FILE',phenolist_filename,.) %>%
              gsub('PHENO_FILE',pheno_filename,.) %>%
              gsub('ANALYSIS_TYPE',tolower(grepl('binary',args$type)),.)
writeLines(json_lines,con=paste0(GWAS_path_tmp,'regenie_',analysis_id,'.json'))

#write phenolist to disk
writeLines(analysis_id,con=paste0(GWAS_path_tmp,phenolist_filename))

#Prepare the phenotype data file
if(args$type=='co_binary'){
  pheno_name <- 'seasonal_val_binary'
}else if(args$type=='co_qt'){
  pheno_name <- 'seasonal_val_qt'
}else if(grepl('cc_binary',args$type)){
  pheno_name <- 'seasonal_val_binary'
}else{
  stop('Type argument not recognized')
}


pheno_dat <- fread('~/Documents/thesis/data/FINNGEN_endpoints_individual_dates.txt')
seasonal_spline <- read_tsv('~/Documents/thesis/data/FINREGISTRY_seasonal_splines.txt')
covariate_dat <- fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz")

seasonal_pheno <- get_seasonal_phenotype(args$endpoint,pheno_dat,seasonal_spline) %>%
                  select(one_of(c('FID',analysis_id=pheno_name)))

if(grepl('co_',args$type)){
  seasonal_pheno <- inner_join(covariate_dat,seasonal_pheno,by='FID')
}else if(grepl('cc_binary_all',args$type)){
  seasonal_pheno <- mutate(seasonal_pheno,seasonal_val_binary=1) %>%
  left_join(covariate_dat,.by='FID') %>%
  mutate(seasonal_val_binary=ifelse(is.na(seasonal_val_binary),0,seasonal_val_binary))
}else if(grepl('high|low',args$type)){
  seasonal_pheno <- left_join(covariate_dat,seasonal_pheno,by='FID') %>%
                    mutate(seasonal_val_binary=ifelse(is.na(seasonal_val_binary),-1,seasonal_val_binary)) %>%
                    filter(if(grepl('high',args$type)) seasonal_val_binary!=0 else seasonal_val_binary!=1) %>%
                    mutate(seasonal_val_binary = as.integer(seasonal_val_binary != -1))
}
seasonal_pheno <- rename(seasonal_pheno,!!analysis_id:=!!pheno_name) %>%
                as.data.table()

fwrite(seasonal_pheno,file=paste0(GWAS_path_tmp,pheno_filename),sep='\t')

# Use gsutils command to copy files to red bucket
base_gsutil_path <- 'gs://fg-production-sandbox-6-red/'

GWAS_path <- paste0(base_gsutil_path,user_path,analysis_id,'/')
system(paste0('gsutil cp -r ',GWAS_path_tmp,' ',GWAS_path))

#Run GWAS pipeline tool
setwd(GWAS_path)
system(paste0("finngen-cli rw -w regenie.wdl -i regenie_",analysis_id,".json -d regenie_sub.zip"))



