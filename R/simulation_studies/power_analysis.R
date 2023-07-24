library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
source('seasonality_dashboard/utils.R')
GWAS_info_dat <- read_tsv('data/finngen_R10_finngen_R10_analysis_data_finngen_R10_pheno_n.tsv') %>%
                 filter(num_gw_significant>1)
seasonal_splines <- read_tsv('seasonality_dashboard/data/FINREGISTRY_seasonal_splines.txt') %>%
                    filter(month == round(month)) %>%
                    mutate(season=ifelse(seasonal_val>0,'high','low'))
monthly_counts_FINNGEN <- read_tsv('data/FINNGEN_endpoints_monthly_count_green.txt') %>%
                          filter_monthly_counts(GWAS_info_dat) %>%
                          inner_join(select(GWAS_info_dat,phenocode,num_controls,num_cases),by=c('ENDPOINT'='phenocode')) %>%
                          inner_join(seasonal_splines,by=c('ENDPOINT','EVENT_MONTH'='month'))

monthly_counts_summary <- group_by(monthly_counts_FINNGEN,ENDPOINT,num_gw_significant,num_cases,num_controls,season) %>%
                          summarise(count=sum(COUNT)) %>%
                          mutate(n=sum(count)) %>%
                          ungroup() %>%
                          pivot_wider(id_cols=c('ENDPOINT','num_cases','num_controls','num_gw_significant','n'),names_from='season',values_from='count')

ggplot(monthly_counts_summary,aes(high,low)) +
geom_point() +
geom_smooth(method='lm') +
geom_abline(intercept = 0,slope=1)


alpha <- 5e-8
thr <- qchisq(alpha,df=1,ncp=0,lower=F)
maf_vec <- seq(0.01,0.5,by=0.01)
beta_vec <- seq(0.01,0.5,by=0.01)
power_dat<- lapply(1:length(beta_vec),function(i){
              lapply(1:length(maf_vec),function(j){
                n_high <- monthly_counts_summary$high
                n_low <- monthly_counts_summary$low
                n <- n_high+n_low
                n_controls <- monthly_counts_summary$num_controls
                n_cases <- monthly_counts_summary$num_cases
                power_vals_cc <-  pchisq(thr,
                                         df = 1,
                                         ncp = (n_high*n_low / n)*2*maf_vec[j]*(1-maf_vec[j])*beta_vec[i]^2, lower = F)
                se_qt <- sqrt((1-2*maf_vec[j]*(1-maf_vec[j])*beta_vec[i]^2)/(2*n*maf_vec[j]*(1-maf_vec[j])))
                power_vals_qt <- pchisq(thr,df=1,ncp=(beta_vec[i]/se_qt)^2,lower=F)
                power_gwas <- pchisq(thr,
                                     df = 1,
                                     ncp = (n_controls*n_cases / (n_cases+n_controls))*2*maf_vec[j]*(1-maf_vec[j])*beta_vec[i]^2, lower = F)
                tibble(ENDPOINT=monthly_counts_summary$ENDPOINT,
                       alpha=alpha,
                       maf=maf_vec[j],
                       beta=beta_vec[i],
                       n_low=n_low,
                       n_high=n_high,
                       n_cases=n_cases,
                       n_controls=n_controls,
                       power_seasonal_cc=power_vals_cc,
                       power_seasonal_qt=power_vals_qt,
                       power_gwas=power_gwas)
              }) %>% bind_rows()
          }) %>% bind_rows()

endpoint_id <- 'F5_DEPRESSIO'
power_by_beta <- filter(power_dat,ENDPOINT==endpoint_id,maf==0.2) %>%
                  pivot_longer(cols = c('power_seasonal_cc','power_seasonal_qt','power_gwas'),names_to='model') %>%
                  mutate(model=gsub('power_','',model)) %>%
                  filter(beta<=0.3) %>%
                  {
                  ggplot(.) +
                  geom_line(aes(beta,value,col=model)) +
                  scale_x_continuous(expand=c(0,0),breaks=seq(0,0.3,by=0.05)) +
                  scale_y_continuous(expand=c(0.01,0)) +
                  xlab(expression(hat(beta))) +
                  ylab('Power') +
                  ggtitle(paste0('MAF = ',unique(.$maf),', n_high = ',unique(.$n_high),', n_low = ',unique(.$n_low))) +
                  theme_bw() +
                  theme(axis.text.x=element_text(size=14),
                        axis.text.y=element_text(size=14),
                        axis.title.x=element_text(size=14),
                        axis.title.y=element_text(size=14),
                        legend.position='none')
                  }

power_by_maf <- filter(power_dat,ENDPOINT==endpoint_id,abs(beta-0.05)<1e-5) %>%
                pivot_longer(cols = c('power_seasonal_cc','power_seasonal_qt','power_gwas'),names_to='model') %>%
                mutate(model=gsub('power_','',model)) %>%
                filter(beta<=0.3) %>%
                {
                  ggplot(.) +
                    geom_line(aes(maf,value,col=model)) +
                    scale_x_continuous(expand=c(0,0),breaks=seq(0,0.5,by=0.1)) +
                    scale_y_continuous(expand=c(0.01,0)) +
                    xlab('MAF') +
                    ylab('') +
                    ggtitle(paste0('beta_hat = ',unique(.$beta),', n_high = ',unique(.$n_high),', n_low = ',unique(.$n_low))) +
                    theme_bw() +
                    theme(axis.text.x=element_text(size=14),
                          axis.text.y=element_text(size=14),
                          axis.title.x=element_text(size=14),
                          axis.title.y=element_text(size=14),
                          legend.text = element_text(size=14),
                          legend.title = element_text(size=16))
                }

power_plot <- power_by_beta + power_by_maf
ggsave('figures/power_heart_failure.png',power_plot,width=12,height=6)

power_by_maf_and_beta <- filter(power_dat,ENDPOINT==endpoint_id) %>%
                          pivot_longer(cols = c('power_seasonal_cc','power_seasonal_qt','power_gwas'),names_to='model') %>%
                          group_by(maf,model) %>%
                          filter(value>=0.9) %>%
                          ungroup() %>%
                          group_by(maf,model) %>%
                          slice(1) %>%
                          mutate(model=gsub('power_','',model)) %>%
                          {
                            ggplot(.) +
                              geom_line(aes(maf,beta,col=model)) +
                              scale_x_continuous(expand=c(0,0),breaks=seq(0,0.5,by=0.1)) +
                              scale_y_continuous(expand=c(0.01,0),limits=c(0,0.3)) +
                              xlab('MAF') +
                              ylab(expression(hat(beta))) +
                              ggtitle(paste0('90% power curves: n_high = ',unique(.$n_high),', n_low = ',unique(.$n_low))) +
                              theme_bw() +
                              theme(axis.text.x=element_text(size=14),
                                    axis.text.y=element_text(size=14),
                                    axis.title.x=element_text(size=14),
                                    axis.title.y=element_text(size=14),
                                    legend.text = element_text(size=14),
                                    legend.title = element_text(size=16))
                          }

ggsave('figures/power_heart_failure_2D.png',power_by_maf_and_beta,width=6,height=4)



