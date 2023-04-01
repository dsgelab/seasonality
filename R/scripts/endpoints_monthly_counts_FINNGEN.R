library(data.table)
library(lubridate)
#Calculate birth date of each individual from detailed longitudinal data. Read first entry of each individual
birth_dat = fread(cmd=paste0("zcat /finngen/library-red/finngen_R10/phenotype_1.0/data/finngen_R10_detailed_longitudinal_1.0.txt.gz ",
                  "| awk -F '\t' 'BEGIN{curr_id=\"\"};{if(curr_id!=$1){curr_id=$1;print $1,$3,$4}}'"))
birth_dates = ymd(birth_dat$APPROX_EVENT_DAY)-dyears(birth_dat$EVENT_AGE)
names(birth_dates) = birth_dat$FINNGENID

path_endpoint = '/finngen/library-red/finngen_R10/phenotype_1.0/data/finngen_R10_endpoint_longitudinal_1.0.txt.gz'
con_endpoint = gzcon(file(path_endpoint,'rb'))
col_names = unlist(strsplit(readLines(con_endpoint,n=1),split='\\t')) #read first line to extract col names
first_event_summary = setDT(NULL)
last_iter_dat = data.table(FINNGENID=character(),ENDPOINT=character())
batch_size = 10e6 # approx 43 iterations in FINNGEN data
time = Sys.time()
while(TRUE){
  dat_endpoint_batch = fread(text=readLines(con_endpoint,n=batch_size),sep='\t')
  if(nrow(dat_endpoint_batch)==0) break #break loop if all data has been read
  names(dat_endpoint_batch) = col_names
  #unique extracts first pair (FINNGENID,ENDPOINT) as table is sorted w.r.t. time
  first_event_batch = unique(dat_endpoint_batch,by=c('FINNGENID','ENDPOINT'))
  event_date = ymd_hms(birth_dates[first_event_batch$FINNGENID])+dyears(first_event_batch$EVENT_AGE)
  first_event_batch$EVENT_MONTH = month(event_date)
  first_event_batch$EVENT_YEAR = year(event_date)
  #anti join with last individual in the previous iteration to avoid double counting (FINNGENID,event) pairs.
  #(necessary due to fixed chunk size)
  first_event_batch = first_event_batch[!last_iter_dt,on=c('FINNGENID','ENDPOINT')]
  #store data on last individual to avoid double counting in next iteration
  last_iter_dat = first_event_batch[FINNGENID==tail(FINNGENID,1),]

  first_event_summary_batch = first_event_batch[,.(COUNT=.N),by=c('ENDPOINT','EVENT_YEAR','EVENT_MONTH')]

  first_event_summary = rbind(first_event_summary,first_event_summary_batch)
  first_event_summary = first_event_summary[,.(COUNT=sum(COUNT)),by=c('ENDPOINT','EVENT_YEAR','EVENT_MONTH')]
}
Sys.time()-time #running time for FINNGEN: 1.17 hours
# Close connection
close(con_endpoint)

# Remove personal data
df_red = first_event_summary
df_green = first_event_summary[first_event_summary$COUNT >= 5, ]

output_path_red = "/data/projects/seasonality/results/FINNGEN_endpoints_monthly_count_red.txt"
output_path_green = "/data/projects/seasonality/results/FINNGEN_endpoints_monthly_count_green.txt"

fwrite(df_red, file = output_path_red, sep = "\t", quote = F)
fwrite(df_green, file = output_path_green, sep = "\t", quote = F)
