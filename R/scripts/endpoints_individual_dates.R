### This script is run within the FinnGen sandbox
library(data.table)
library(lubridate)
freeze <- 'R11'
detailed_pheno_path <- paste0('/finngen/library-red/finngen_',freeze,'/phenotype_1.0/data/finngen_',freeze,'_detailed_longitudinal_1.0.txt.gz')
#Calculate birth date of each individual from detailed longitudinal data. Read first entry of each individual
birth_dat <- fread(cmd=paste0("zcat ",detailed_pheno_path," ",
                             "| awk -F '\t' 'BEGIN{curr_id=\"\"};{if(curr_id!=$1){curr_id=$1;print $1,$3,$4}}'"))
birth_dates <- ymd(birth_dat$APPROX_EVENT_DAY)-dyears(birth_dat$EVENT_AGE)
names(birth_dates) <- birth_dat$FINNGENID
GWAS_info_dat <- fread(paste0('/finngen/library-green/finngen_',freeze,'/finngen_',freeze,'_analysis_data/finngen_',freeze,'_pheno_n.tsv'))
#restrict to case-control phenotypes with sufficiently many cases and some significant GWAS signal
GWAS_info_dat <- GWAS_info_dat[num_gw_significant>0 & num_controls>0 & num_cases>5000,'phenocode']
path_endpoint = paste0('/finngen/library-red/finngen_',freeze,'/phenotype_1.0/data/finngen_',freeze,'_endpoint_longitudinal_1.0.txt.gz')
con_endpoint = gzcon(file(path_endpoint,'rb'))
col_names = unlist(strsplit(readLines(con_endpoint,n=1),split='\\t')) #read first line to extract col names
first_event_dat = setDT(NULL)
batch_size = 10e6 # approx 43 iterations in FINNGEN data
time = Sys.time()
i <- 1
while(TRUE){
  print(i)
  dat_endpoint_batch = fread(text=readLines(con_endpoint,n=batch_size),sep='\t')
  if(nrow(dat_endpoint_batch)==0) break #break loop if all data has been read
  names(dat_endpoint_batch) = col_names
  dat_endpoint_batch <- dat_endpoint_batch[GWAS_info_dat,on=c('ENDPOINT'='phenocode'), nomatch=0]
  #unique extracts first pair (FINNGENID,ENDPOINT) as table is sorted w.r.t. time
  first_event_dat_batch <- unique(dat_endpoint_batch,by=c('FINNGENID','ENDPOINT'))
  first_event_dat_batch$EVENT_DATE <- as_date(ymd_hms(birth_dates[first_event_dat_batch$FINNGENID]) +
                                      dyears(first_event_dat_batch$EVENT_AGE))
  first_event_dat_batch$EVENT_YEAR <- year(first_event_dat_batch$EVENT_DATE)
  first_event_dat <- rbind(first_event_dat,first_event_dat_batch)
  i <- i + 1
}
#Extract first element per (FINNGENID,ENDPOINT) pair, in case there are duplicates due to batching
first_event_dat <- unique(dat_endpoint_batch,by=c('FINNGENID','ENDPOINT'))
Sys.time()-time #running time for FINNGEN: 1.17 hours
# Close connection
close(con_endpoint)

output_path = paste0("/data/FINNGEN_",freeze,"_endpoints_individual_dates.txt")
fwrite(first_event_dat, file = output_path, sep = "\t", quote = F)
