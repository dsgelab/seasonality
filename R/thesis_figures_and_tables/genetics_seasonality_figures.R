library(qqman)
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(xtable)

#### Help functions ####
summarise_top_signals <- function(e,type){
  path <- paste0('results/GWAS_results/',e,'.',type,'_summary.txt')
  most_severe_lab <- c('intergenic_variant'='Intergenic variant','intron_variant'='Intron variant',
                       '3_prime_UTR_variant'="3' UTR variant",'missense_variant'='Missense variant',
                       'upstream_gene_variant'='Upstream gene variant')
  top_signal_summary <- read_tsv(path) %>%
    dplyr::filter(pval<5e-8) %>%
    group_by(`#chrom`) %>%
    slice(which.max(mlogp)) %>%
    arrange(pval) %>%
    ungroup() %>%
    mutate(endpoint=e,
           OR=exp(beta),
           gene_most_severe=as.character(ifelse(is.na(gene_most_severe),'',gene_most_severe)),
           most_severe=most_severe_lab[most_severe],
           rsid=paste0(rsid,'[',alt,']')) %>%
    select(endpoint,`#chrom`,rsid,pval,beta,OR,af_alt,gene_most_severe,most_severe)
}

plot_manhattan <- function(file,snps,endpoint_name,pval_thres=0.001,flank=1e5,...){
  sumstats <- fread(file)
  sumstats_signif <- sumstats[pval<pval_thres,]
  #random_idx <- sample(which(sumstats$pval>pval_thres),size = 20000)
  #sumstats_random <- sumstats[random_idx,]
  #sumstats_signif <- rbind(sumstats_signif,sumstats_random)
  sumstats_signif$SNP <- 1:nrow(sumstats_signif)
  snps_of_interest <- rep(F,nrow(sumstats_signif))
  for(snp in snps){
    sumstats_signif$SNP <- ifelse(sumstats_signif$`#chrom`== snp$chrom & sumstats_signif$pos==snp$pos,snp$name,sumstats_signif$SNP)
    snps_of_interest <- snps_of_interest | (sumstats_signif$`#chrom`==snp$chrom & sumstats_signif$pos > snp$pos-flank & sumstats_signif$pos < snp$pos+flank)
  }
  manhattan(sumstats_signif,chr='#chrom',bp='pos',p='pval',highlight=sumstats_signif$SNP[snps_of_interest],suggestiveline=FALSE,annotatePval = 5e-8,chrlabs=c(1:22,'X'),...)
  title(main=endpoint_name,line=1)
}

#### List of endpoints and the number of cases in summer and winter ####
endpoint_names <- read_excel('R/thesis_figures_and_tables/endpoint_name_mapping.xlsx')

phenotype_summary <- bind_rows(read_tsv('results/GWAS_results/phenotype_summary.txt'),
                               read_tsv('results/GWAS_results/phenotype_summary_infectious.txt'),
                               read_tsv('results/GWAS_results/phenotype_summary_auto_infl.txt')) %>%
  separate(analysis_id,into=c('endpoint','mod_type'),sep='\\.') %>%
  filter(mod_type=='co_binary') %>%
  inner_join(endpoint_names,.,by=c('Endpoint ID'='endpoint')) %>%
  arrange(desc(n)) %>%
  mutate_at(c('n','cases','controls'),function(x) ifelse(nchar(x)>3,paste0(substr(x,1,nchar(x)-3),',',substr(x,nchar(x)-2,nchar(x))),x))
print(xtable(select(phenotype_summary,Endpoint,n,`n high season`=cases,`n low season`=controls),digits=c(0,0,0,0,0)),include.rownames=FALSE)

#### Association summary for genome-wide significant variants ######
top_signals_summary <- lapply(phenotype_summary$`Endpoint ID`,function(e) {
  summarise_top_signals(e=e,type='co_binary')
}) %>% bind_rows() %>%
  inner_join(endpoint_names,.,by=c('Endpoint ID'='endpoint')) %>%
  arrange(pval) %>%
  mutate(pval=formatC(pval,format='e',digits=1),
         OR=format(OR,digits=2,nsmall=2),
         af_alt=paste0(format(round(100*af_alt,1),nsmall=1),'%'),
         gene_most_severe=ifelse(nchar(gene_most_severe)>0,paste0('\\textit{',gene_most_severe,'}'),gene_most_severe))


print(xtable(select(top_signals_summary,Endpoint,rsID=rsid,`MAF`=af_alt,Gene=gene_most_severe,`P-value`=pval,OR)),include.rownames=FALSE)

##### Manhattan plots ######
filename <- 'pdf/manhattan_plots.pdf'
pdf(file = filename, width = 7, height = 7)
par(mfrow=c(3,1),oma=c(0.1,1.1,0.1,3.1))
plot_manhattan(file='results/GWAS_results/F5_DEPRESSIO.co_binary.gz',snps=list(list(name='rs149442832',chrom=6,pos=76678973),
                                                                               list(name='rs75529310',chrom=7,pos=20412850)),
               endpoint_name='Major depression',pval_thres=0.001,flank=1e5,xlab='',mar=c(0.1,4.1,1.1,2.1),ylim=c(3,8.5))
plot_manhattan(file='results/GWAS_results/AB1_INTESTINAL_INFECTIONS.co_binary.gz',snps=list(list(name='rs516246',chrom=19,pos=48702915)),
               endpoint_name='Intestinal infections',pval_thres=0.001,flank=1e5,xlab='',mar=c(0.1,4.1,1.1,2.1),ylim=c(3,8.5))
plot_manhattan(file='results/GWAS_results/C_STROKE.co_binary.gz',snps=list(list(name='rs1408687730',chrom=3,pos=66766028)),
               endpoint_name='Stroke',pval_thres=0.001,flank=1e5,mar=c(5.1,4.1,1.1,2.1),ylim=c(3,8.5))
dev.off()




##### Sensitivity analysis #####
monthly_counts_FinRegistry  <- read_tsv('seasonality_dashboard/data/FINREGISTRY_endpoints_monthly_count_green.txt') %>%
  filter_monthly_counts(.,GWAS_info_dat=GWAS_info_dat,endpoint_mapping_path='seasonality_dashboard/data/finngen_R10_endpoint_core_noncore_1.0.xlsx')
endpoints=c('F5_DEPRESSIO','AB1_INTESTINAL_INFECTIONS','C_STROKE')
a <- 1;b <- 13
sensitivity_spline_dat <- lapply(endpoints, function(e){
                          seasonal_splines_endpoint <- extract_seasonal_spline(filter(monthly_counts_FinRegistry,ENDPOINT==e),a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=6)
                          seasonal_splines_endpoint_k4 <- extract_seasonal_spline(filter(monthly_counts_FinRegistry,ENDPOINT==e),a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=4)
                          seasonal_splines_endpoint_k8 <- extract_seasonal_spline(filter(monthly_counts_FinRegistry,ENDPOINT==e),a=a,b=b,adjustment='',seasonal_spline_type='cp',k_seasonal=8)
                          seasonal_splines_endpoint_adj <- extract_seasonal_spline(filter(monthly_counts_FinRegistry,ENDPOINT==e),a=a,b=b,adjustment='binary',seasonal_spline_type='cp',k_seasonal=6)
                          bind_rows(list(k6=seasonal_splines_endpoint,
                                         k4=seasonal_splines_endpoint_k4,
                                         k8=seasonal_splines_endpoint_k8,
                                         adjusted=seasonal_splines_endpoint_adj),.id='spline_type') %>%
                          mutate(ENDPOINT=e) %>% select(ENDPOINT,everything())
                        }) %>% bind_rows() %>%
                        mutate(ENDPOINT=factor(ENDPOINT,levels=endpoints,labels=c('Major depression','Intestinal infections','Stroke')),
                               spline_type=factor(spline_type,levels=c('k6','k4','k8','adjusted'),labels=c('Used','Knots=4','Knots=8','Adjustment')))


sensitivity_plot <- ggplot(sensitivity_spline_dat,aes(month,seasonal_val,col=spline_type)) +
                    geom_line(linewidth=1.2) +
                    facet_wrap(~ENDPOINT,scales='free') +
                    geom_hline(yintercept=0,linetype='dashed') +
                    scale_color_manual(values=c("#000000","#BC3C29FF", "#0072B5FF", "#E18727FF"),name='Analysis type') +
                    scale_x_continuous(breaks=seq(1,12,by=2)) +
                    xlab('Month number') +
                    ylab('Smooth term') +
                    theme_classic() +
                    theme(strip.text = element_text(size = 12, margin = margin(),hjust = 0),
                          strip.background = element_blank())
save_tikz_plot(sensitivity_plot,width=7,height=3,filename = 'tex/sensitivity_plot.tex')

sensitivity_table <- read_excel('results/sensitivity_analysis.xlsx') %>%
                     mutate(OR=exp(beta)) %>%
                     mutate(pval=formatC(pval,format='e',digits=2),
                            OR=format(OR,digits=2,nsmall=2)) %>%
                    select(Endpoint,`Analysis type`,OR,`p-value`=pval)
print(xtable(sensitivity_table),include.rownames=FALSE)


##### Supplementary table with the mapping between endpoint ids and #####
endpoint_name_mapping <- read_excel('R/thesis_content/endpoint_name_mapping.xlsx')
print(xtable(endpoint_name_mapping),include.rownames=F)
