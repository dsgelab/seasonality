library(readr)
library(dplyr)
library(tidyr)

var_compare <- c('type','adjustment')
freezes <- c('R10','R11')
base_results_path <- 'seasonality_paper/results/'

job_info <- lapply(freezes,function(freeze) read_tsv(paste0(base_results_path,'GWAS_',freeze,'/','job_info.tsv')))
names(job_info) <- freezes
results_summary <- lapply(freezes,function(freeze){
                      analysis_ids <- job_info[[freeze]]$analysis_id
                      lapply(analysis_ids,function(a_id){
                        read_tsv(paste0(base_results_path,'GWAS_',freeze,'/',a_id,'/',a_id,'_summary.txt'),show_col_types=F) %>%
                        select(-fin.AF,-nfsee.AF) %>%
                        mutate(analysis_id=a_id)
                      }) %>% bind_rows() %>% mutate(freeze=freeze)
                    }) %>% bind_rows() %>%
                    mutate(ENDPOINT=gsub('_co_.*','',analysis_id),
                           analysis_id_suffix=sapply(1:n(), function(i) gsub(paste0(ENDPOINT[i],'_'),'',analysis_id[i]))) %>%
                    separate(analysis_id_suffix,sep='_',into=c('design','type','adjustment','time_period')) %>%
                    select(ENDPOINT,type,adjustment,everything())

signif_remaining <- arrange(results_summary,pval)
top_table <- tibble()
while(nrow(signif_remaining)>0){
  top_variant <- slice(signif_remaining,1)
  signif_remaining <- filter(signif_remaining,!(ENDPOINT==top_variant$ENDPOINT &
                                                type==top_variant$type &
                                                adjustment==top_variant$adjustment &
                                                pos < top_variant$pos + 5e5 &
                                                pos > top_variant$pos - 5e5))
  top_table <- bind_rows(top_table,top_variant)
}


analysis_names <- c('co_binary_unadj','co_binary_adj','co_qt_unadj','co_qt_adj')
tmp_dir <- tempdir()
GWAS_results_other <- lapply(freezes,function(fr){
                        lapply(analysis_ids,function(a_id){
                          sumstat_path <- job_info[[fr]]$sumstat_path[job_info[[fr]]$analysis_id==a_id]
                          sumstat_index_path <- job_info[[fr]]$sumstat_index_path[job_info[[fr]]$analysis_id==a_id]
                          file.copy(sumstat_path,tmp_dir)
                          file.copy(sumstat_index_path,tmp_dir)
                          endpoint <- gsub('_co_.*','',a_id)
                          top_table_others <- filter(top_table,analysis_id!=a_id | freeze==fr,ENDPOINT==endpoint) %>%
                                                    distinct(`#chrom`,pos,REF,ALT,rsid)
                          results_other_pheno <- lapply(1:nrow(top_table_others),function(i){
                                                  results_raw <- system(paste0('tabix ',tmp_dir,'/',basename(sumstat_path),' ',
                                                                               top_table_others$`chrom`[i],':',
                                                                               top_table_others$pos[i],'-',top_table_other  s$pos[i]),intern=T)
                                                  if(!identical(results_raw,character(0))){
                                                    read.table(text=results_raw,fill=T) %>%
                                                    as_tibble() %>%
                                                    mutate(V3=if(is.logical(V3)) 'T' else V3,
                                                           V4=if(is.logical(V4)) 'T' else V4) %>%
                                                    select(1:9) %>%
                                                    filter(V3==top_table_others$ref[i],V4==top_table_others$alt[i])
                                                  }else{
                                                    tibble(V1=as.integer(top_table_others$`#chrom`[i]),V2=as.integer(top_table_others$pos[i]),
                                                           V3=top_table_others$ref[i],V4=top_table_others$alt[i],
                                                           V5=as.numeric(NA),V6=as.numeric(NA),
                                                           V7=as.numeric(NA),V8=as.numeric(NA),
                                                           V9=as.numeric(NA))
                                                  }
                                            }) %>% bind_rows() %>%
                                            mutate(analysis_id=a_id)
                          names(results_other_pheno) <- c('#chrom','pos','ref','alt','pval','mlogp','beta','sebeta','af_alt','analysis_id')
                          file.remove(paste0(tmp_dir,'/',basename(sumstats_path)))
                          file.remove(paste0(tmp_dir,'/',basename(sumstats_index_path)))
                          return(results_other_pheno)
                        }) %>% bind_rows() %>% mutate(freeze=fr)
                      }) %>% bind_rows()

cmp_dat <- bind_rows(select(top_table,ENDPOINT,type,adjustment,feeze,`#chrom`,pos,ref,alt,pval,mlogp,beta,sebeta,af_alt),
                     select(GWAS_results_other,ENDPOINT,type,adjustment,feeze,`#chrom`,pos,ref,alt,pval,mlogp,beta,sebeta,af_alt)) %>%
           distinct(ENDPOINT,type,adjustment,feeze,`#chrom`,pos,ref,alt,.keep_all=T)

filter(cmp_dat,af_alt>0.01,adjustment=='unadj',freeze=='R11') %>%
pivot_wider(id_cols=c('ENDPOINT','#chrom','pos','ref','alt'),names_from='type',values_from=c('mlogp')) %>%
ggplot(aes(x=)) +
geom_point()
