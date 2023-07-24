library(dplyr)
library(readr)
library(lubridate)
library(data.table)
source('utils.R')

prepare_seasonality_phenotype <- function(endpoint,type){
  #Create folder structure for GWAS
  base_path <- '/finngen/red/'
  user_path <- 'solviro/GWAS_R11/'
  analysis_id <- paste(endpoint,type,sep='_')
  GWAS_path_tmp <- paste0('~/Documents/thesis/tmp/')
  template_path <- paste0(base_path,user_path,'template/')
  phenolist_filename <-  paste0('pheno_',analysis_id,'.txt')
  pheno_filename <-  paste0(analysis_id,'.txt.gz')
  #create directory of not yet created
  if(!dir.exists(GWAS_path_tmp)){
    dir.create(GWAS_path_tmp)
  }

  #write phenolist to disk
  writeLines(analysis_id,con=paste0(GWAS_path_tmp,phenolist_filename))

  #Prepare the phenotype data file
  if(type=='co_binary'){
    pheno_name <- 'seasonal_val_binary'
  }else if(type=='co_qt'){
    pheno_name <- 'seasonal_val_qt'
  }else if(grepl('cc_binary',type)){
    pheno_name <- 'seasonal_val_binary'
  }else{
    stop('Type argument not recognized')
  }


  pheno_dat <- fread('~/Documents/thesis/data/FINNGEN_endpoints_individual_dates.txt')
  seasonal_spline <- read_tsv('~/Documents/thesis/data/FINREGISTRY_seasonal_splines.txt')
  covariate_dat <- fread("/finngen/library-red/finngen_R11/analysis_covariates/R11_COV_V0.FID.txt.gz")

  adj_cov = strsplit(paste0("SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,",
                            "BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,",
                            "BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm_relift,BATCH_DS7_FINRISK_FR12_norm,",
                            "BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm,BATCH_DS10_FINRISK_Palotie_norm,",
                            "BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,",
                            "BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,",
                            "BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm_relift,BATCH_DS18_MIGRAINE_1_norm_relift,",
                            "BATCH_DS19_MIGRAINE_2_norm,BATCH_DS20_SUPER_1_norm_relift,BATCH_DS21_SUPER_2_norm_relift,",
                            "BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm_nosymmetric,BATCH_DS24_SUPER_3_norm,",
                            "BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS26_DIREVA_norm,BATCH_DS27_NFBC66_norm,BATCH_DS28_NFBC86_norm"),
                     split=',') %>% unlist() %>% c(.,paste0('PC',1:10))

  seasonal_pheno <- get_seasonal_phenotype(endpoint,pheno_dat,seasonal_spline) %>%
    select(one_of(c('FID',analysis_id=pheno_name)))

  if(grepl('co_',type)){
    seasonal_pheno <- inner_join(covariate_dat,seasonal_pheno,by='FID')
  }else if(grepl('cc_binary_all',type)){
    seasonal_pheno <- mutate(seasonal_pheno,seasonal_val_binary=1) %>%
      left_join(covariate_dat,.by='FID') %>%
      mutate(seasonal_val_binary=ifelse(is.na(seasonal_val_binary),0,seasonal_val_binary))
  }else if(grepl('high|low',type)){
    seasonal_pheno <- left_join(covariate_dat,seasonal_pheno,by='FID') %>%
      mutate(seasonal_val_binary=ifelse(is.na(seasonal_val_binary),-1,seasonal_val_binary)) %>%
      filter(if(grepl('high',type)) seasonal_val_binary!=0 else seasonal_val_binary!=1) %>%
      mutate(seasonal_val_binary = as.integer(seasonal_val_binary != -1))
  }
  seasonal_pheno <- rename(seasonal_pheno,!!analysis_id:=!!pheno_name) %>%
    as.data.table()

  #check if any batch has one or fewer individuals -> this will lead to REGENIE failure
  batch_dat=select(seasonal_pheno,matches(adj_cov[grepl('BATCH',adj_cov)]))
  small_batches = names(batch_dat)[colSums(batch_dat)<=1]
  if(length(small_batches) != 0){
    warning(paste0('The following batches are too small to be included in the GWAS:\n',
                   paste(small_batches,collapse='\n')))
  }

  #Edit json file inline with the current GWAS to be run
  json_lines <- readLines(paste0(template_path,'regenie_ENDPOINT.json')) %>%
    gsub('USER_PATH',paste0(user_path,analysis_id,'/'),.) %>%
    gsub('PHENOLIST_FILE',phenolist_filename,.) %>%
    gsub('PHENO_FILE',pheno_filename,.) %>%
    gsub('ANALYSIS_TYPE',tolower(grepl('binary',type)),.) %>%
    gsub('ADJUSTMENT_COVARIATES',paste(adj_cov[!(adj_cov %in% small_batches)],collapse=','),.)
  writeLines(json_lines,con=paste0(GWAS_path_tmp,'regenie_',analysis_id,'.json'))

  fwrite(seasonal_pheno,file=paste0(GWAS_path_tmp,pheno_filename),sep='\t')

  # Use gsutils command to copy files to red bucket
  base_gsutil_path <- 'gs://fg-production-sandbox-6-red/'
  GWAS_gsutil_path <- paste0(base_gsutil_path,user_path,analysis_id,'/')
  system(paste0('gsutil cp ',GWAS_path_tmp,'/* ',GWAS_gsutil_path))

  unlink(GWAS_path_tmp,recursive=T)
}


#Run GWAS pipeline tool
# GWAS_local_path <- paste0('/finngen/red/',user_path,analysis_id,'/')
# setwd(GWAS_local_path)
# system(paste0("finngen-cli rw -w regenie.wdl -i regenie_",analysis_id,".json -d regenie_sub.zip"))



