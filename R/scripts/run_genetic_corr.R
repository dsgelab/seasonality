#!/usr/bin/env Rscript
library(data.table)
library(optparse)

munge_sumstats <- function(sumstats_path,HM3_dat,ldscref,ldsc_build,suffix){
  print('Reading sumstats file')
  sumstats <- fread(sumstats_path,select=c('#chrom','pos','ref','alt','pval','beta','sebeta'),
                     col.names=c('CHR','BP','REF','ALT','P','BETA','SE'))
  sumstats[,CHR := paste0('chr',CHR)]
  sumstats_hm3_match <- merge(sumstats,HM3_dat,by=c('CHR','BP','REF','ALT'))
  sumstats_hm3_antimatch <- merge(sumstats,HM3_dat,by.x=c('CHR','BP','REF','ALT'),by.y=c('CHR','BP','REF','ALT'))
  sumstats_hm3 <- rbind(sumstats_hm3_match,sumstats_hm3_antimatch)
  n <- nrow(sumstats_hm3)
  tmp_dir <- tempdir()
  print(paste0('Writing into tmpdir: ',tmp_dir))
  sumstats_path <- paste0(tmp_dir,'/sumstats_tmp',suffix,'.tsv')
  fwrite(sumstats_hm3[,hg38:=NULL],file = sumstats_path,sep='\t')
  munged_outpath <- paste0(tmp_dir,'/munged',suffix)
  system(paste0(ldsc_build,'munge_sumstats.py',
                    ' --sumstats ',sumstats_path,
                    ' --out ',munged_outpath,
                    ' --N ',n,
                    ' --a1 ALT',
                    ' --a2 REF',
                    ' --p P',
                    ' --signed-sumstats BETA,0',
                    ' --merge-alleles ',ldscref,'/w_hm3.snplist',
                    ' --chunksize 500000'))
  return(paste0(munged_outpath,'.sumstats.gz'))
}

run_genetic_correlation <- function(sumstats1,sumstats2,outpath,hm3,ldscref,ldsc_build){
  print('Reading HM3 file:')
  HM3_dat <- fread(hm3,select=c('SNP','hg38','REF','ALT'))
  HM3_dat[, c("CHR", "BP") := tstrsplit(hg38, ":", fixed=TRUE)]
  HM3_dat[,BP:=as.integer(BP)]

  print('Munging ongoing')
  munged_sumstats1_path <- munge_sumstats(sumstats_path=sumstats1,HM3_dat=HM3_dat,ldscref=ldscref,ldsc_build=ldsc_build,suffix='1')
  munged_sumstats2_path <- munge_sumstats(sumstats_path=sumstats2,HM3_dat=HM3_dat,ldscref=ldscref,ldsc_build=ldsc_build,suffix='2')

  system(paste0(ldsc_build,'ldsc.py',
                ' --rg ',munged_sumstats1_path,',',munged_sumstats2_path,
                ' --ref-ld-chr ',ldscref,'/fin_ldsc/',
                ' --w-ld-chr ',ldscref,'/fin_ldsc/',
                ' --out ',outpath))
  unlink(basename(dirname(munged_sumstats1_path))) # remove temp directory
}

run_seasonal_genetic_correlation <- function(endpoint,type){
  hm3 <- '~/Documents/thesis/ldsc_calc/HM3Ref'
  ldscref <- '/finngen/library-green/ldsc'
  ldsc_build <- '~/Documents/thesis/ldsc_calc/ldsc/'
  outfolder <- '~/Documents/thesis/GWAS_R10/results/genetic_correlation/'

  if(type=='binary_vs_qt'){
    sumstats1 <- paste0('~/Documents/thesis/GWAS_R10/results/',endpoint,'.co_binary/binary.gz')
    sumstats2 <- paste0('~/Documents/thesis/GWAS_R10/results/',endpoint,'.co_qt/qt.gz')
    outname <- paste0(endpoint,'.binary_vs_qt.corr')
  }else if(type=='seasonal_vs_cc'){
    sumstats1 <- paste0('~/Documents/thesis/GWAS_R10/results/',endpoint,'.co_binary/binary.gz')
    sumstats2 <- paste0('/finngen/library-green/',endpoint,'.gz')
    outname <- paste0(endpoint,'.seasonal_vs_cc.corr')
  }
  outpath <- paste0(outfolder,outname)
  run_genetic_correlation(sumstats1=sumstats1,sumstats2=sumstats2,outpath=outpath,hm3=hm3,ldscref=ldscref,ldsc_build=ldsc_build)
}

endpoints <- c('AB1_INTESTINAL_INFECTIONS','AB1_GASTROENTERITIS_NOS',
               'C_STROKE','F5_DEPRESSIO','H7_CONJUNCTIVITIS','H8_MIDDLEMASTOID',
               'J10_INFLUPNEU','J10_PNEUMONIA','J10_SINUSITIS')
for(e in endpoints){
  run_seasonal_genetic_correlation(e,type='binary_vs_qt')
}



