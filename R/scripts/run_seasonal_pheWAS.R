library(dplyr)
library(tidyr)
get_seasonal_phenotype <- function(endpoint_ids,pheno_dat,seasonal_spline,spline_endpoint=NULL){
  dates <- seq(ymd('1998-01-01'),ymd('2019-12-31'),by='1 day')
  dates_dat <- tibble(year=year(dates),month=month(dates),day=day(dates)) %>%
                group_by(year,month) %>%
                summarise(num_days=max(day))
  pheno_dat_endpoint <- as_tibble(pheno_dat[ENDPOINT %in% endpoint_ids,]) %>%
                        rename(FID=FINNGENID) %>%
                        mutate(EVENT_MONTH=month(EVENT_DATE),
                               EVENT_DAY=day(EVENT_DATE)) %>%
                        inner_join(dates_dat,by=c('EVENT_YEAR'='year','EVENT_MONTH'='month')) %>%
                        mutate(EVENT_MONTH_DEC=EVENT_MONTH + EVENT_DAY/num_days) %>%
                        mutate(EVENT_MONTH_DEC=ifelse(EVENT_MONTH_DEC>12.5,
                                                      EVENT_MONTH_DEC-12,
                                                      EVENT_MONTH_DEC))
  if(is.null(splice_endpoint)){

  }
  seasonal_spline_endpoint <- filter(seasonal_spline,ENDPOINT==endpoint_id)
  pheno_dat_endpoint <- mutate(pheno_dat_endpoint,
                               seasonal_val=approx(x=seasonal_spline_endpoint$month,
                                                   y=seasonal_spline_endpoint$seasonal_val,
                                                   xout=EVENT_MONTH_DEC)$y,
                               seasonal_val_qt=qnorm(rank(seasonal_val)/(n()+0.5)),
                               seasonal_val_binary=as.integer(seasonal_val > 0))
  return(pheno_dat_endpoint)
}



chrom <- 'chr7'
pos <- 20412850
gt_dat <- system(paste0("tabix -h /finngen/library-red/finngen_R10/genotype_1.0/data/finngen_R10_",
                        chrom,".vcf.gz ",chrom,":",pos,"-",pos,
                        " | awk '$1 !~ \"##\"' | datamash transpose | awk '$1 ~ \"^FG\"'"),intern=T) %>%
          read.table(text=.,header=F) %>%
          separate(V2,sep=':',into=c('GT','DS','GP')) %>%
          separate(GT,sep='\\|',into=c('h1','h2')) %>%
          mutate(h1=as.numeric(h1),
                 h2=as.numeric(h2),
                 gt=h1+h2) %>%
          rename(FID=V1)

pheno_dat <- fread('~/Documents/thesis/data/FINNGEN_endpoints_individual_dates.txt')
seasonal_spline <- read_tsv('~/Documents/thesis/data/FINREGISTRY_seasonal_splines.txt')
covariate_dat <- fread("/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz")

seasonal_pheno <- get_seasonal_phenotype(args$endpoint,pheno_dat,seasonal_spline) %>%
  select(one_of(c('FID',analysis_id=pheno_name)))
