library(dplyr)
library(tidyr)
library(data.table)
library(lubridate)
source('utils.R')
chrom <- 'chr7'
pos <- 20412850
spline_endpoint='F5_DEPRESSIO'
spline_type='k4'
#chrom= 'chr19'
#pos=48703417
#spline_endpoint='AB1_INTESTINAL_INFECTIONS'
#chrom= 'chr3'
#pos=66766028
#spline_endpoint='C_STROKE'

# extract genotype for the input variant
gt_dat <- get_gt_dat(chrom=chrom,pos=pos)

# read in necessary files
pheno_dat <- fread('~/Documents/thesis/data/FINNGEN_endpoints_individual_dates.txt')
seasonal_spline <- fread(paste0('~/Documents/thesis/data/FINREGISTRY_seasonal_splines_',spline_type,'.txt'))
covariate_dat <- fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz")

#R10 adjustment covariates
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

# Run association to all seasonal phenotypes
pheWAS_result_table = lapply(unique(seasonal_spline$ENDPOINT),function(e){
                        seasonal_pheno <- get_seasonal_phenotype(e,pheno_dat,seasonal_spline,spline_endpoint = if(!is.null(spline_endpoint)) spline_endpoint else e)
                        if(nrow(seasonal_pheno)>0){
                          full_dat=inner_join(covariate_dat,seasonal_pheno,by='FID') %>%
                            inner_join(gt_dat,by='FID') %>%
                            select(one_of(c('FID','seasonal_val_binary','seasonal_val_01','seasonal_val_qt','gt',adj_cov)))
                          f_base=paste0('gt + ',paste(adj_cov[adj_cov %in% names(full_dat)],collapse='+'))
                          f_binary=as.formula(paste0('seasonal_val_binary ~ ',f_base))
                          mod_binary=glm(f_binary,data=full_dat,family='binomial')
                          summary_binary=summary(mod_binary)$coefficients[2,c(1,2,4)] %>% 
                                         as.list() %>%
                                         as_tibble()
                          names(summary_binary)=c('est_binary','se_binary','pval_binary')
                          f_continuous=as.formula(paste0('seasonal_val_01 ~ ',f_base))
                          mod_continuous=glm(f_continuous,data=full_dat,family='gaussian')
                          summary_continuous=summary(mod_continuous)$coefficients[2,c(1,2,4)] %>% 
                                              as.list() %>%
                                              as_tibble()
                          names(summary_continuous)=c('est_continuous','se_continuous','pval_continuous')
                          f_qt=as.formula(paste0('seasonal_val_qt ~ ',f_base))
                          mod_qt=glm(f_qt,data=full_dat,family='gaussian')
                          summary_qt=summary(mod_qt)$coefficients[2,c(1,2,4)] %>% 
                                      as.list() %>%
                                      as_tibble()
                          names(summary_qt)=c('est_qt','se_qt','pval_qt')
                          bind_cols(summary_binary,summary_continuous,summary_qt) %>%
                          mutate(endpoint=e) %>%
                          select(endpoint,everything())
                        }else{
                          NULL
                        }
                      
                      }) %>% bind_rows()
write.table(pheWAS_result_table,file=paste0('~/Documents/thesis/pheWAS/',chrom,'_',pos,'_',ifelse(is.null(spline_endpoint),'no',''),'match.txt'),row.names=F,quote=F,sep='\t')
