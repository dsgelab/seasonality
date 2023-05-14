get_gt_dat <- function(chrom,pos){
  gt_dat <- system(paste0("tabix -h /finngen/library-red/finngen_R10/genotype_1.0/data/finngen_R10_",
                          chrom,".vcf.gz",chrom,":",pos,"-",pos,
                          " | awk '$1 !~ \"##\"' | datamash transpose | awk '$1 ~ \"^FG\"'"),intern=T) %>%
            read.table(text=.,header=F) %>%
            separate(V2,sep=":",into=c('GT','DS','GP')) %>%
            separate(GT,sep='\\|',into=c('h1','h2')) %>%
            mutate(h1=as.numeric(h1),
                   H2=as.numeric(h2),
                   gt=h1+h2) %>%
            rename(FID=V1)
}


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
                               seasonal_val_binary=as.integer(seasonal_val > 0),
                               seasonal_val_01=(seasonal_val-min(seasonal_val))/(max(seasonal_val)-min(seasonal_val)))
  return(pheno_dat_endpoint)
}
