get_endpoint_mapping <- function(endpoint_vec,endpoint_categories){
  endpoint_groups <- mutate(endpoint_categories,code=gsub('#','',code)) %>%
    pivot_longer(cols=c('CHAPTER','OTHER'),names_to='descr_type',values_to='descr',values_drop_na = T)
  endpoint_mapping <- mutate(tibble(ENDPOINT=endpoint_vec),prefix=gsub('_.*$','',ENDPOINT)) %>%
                      select(ENDPOINT,prefix) %>%
                      left_join(endpoint_groups,by=c('prefix'='code'))
  return(endpoint_mapping)
}

filter_monthly_counts <- function(monthly_counts,year_start,year_end){
  monthly_counts_filtered <- filter(monthly_counts,EVENT_YEAR>=year_start & EVENT_YEAR<=year_end)
  month_year_endpoint_template <- expand_grid(ENDPOINT=unique(monthly_counts_filtered$ENDPOINT),
                                              EVENT_YEAR=seq(min(monthly_counts_filtered$EVENT_YEAR),
                                                             max(monthly_counts_filtered$EVENT_YEAR)),
                                              EVENT_MONTH=seq_len(12))
  #endpoint_mapping <- get_endpoint_mapping(path=endpoint_mapping_path,seasonal_summary=distinct(monthly_counts_filtered,ENDPOINT))
  monthly_counts_filtered <- left_join(month_year_endpoint_template,
                                       monthly_counts_filtered,
                                       by=c('ENDPOINT','EVENT_YEAR','EVENT_MONTH')) %>%
    mutate(COUNT=ifelse(is.na(COUNT),0L,COUNT),
           EVENT_DATE=ym(paste0(EVENT_YEAR,'-',EVENT_MONTH))) %>%
    #left_join(select(endpoint_mapping,ENDPOINT,prefix,brief),by='ENDPOINT') %>%
    #filter(is.na(prefix) | !(prefix %in% c('OTHER','VWXY20'))) %>%
    arrange(ENDPOINT,EVENT_YEAR,EVENT_MONTH)
  #all endpoints have to have more than 0 diagnoses each year.
  monthly_counts_filtered  <- group_by(monthly_counts_filtered,ENDPOINT,EVENT_YEAR) %>%
    summarise(count=sum(COUNT)) %>%
    group_by(ENDPOINT) %>%
    filter(all(count>0)) %>%
    distinct(ENDPOINT) %>%
    inner_join(monthly_counts_filtered,.,by='ENDPOINT')
  return(monthly_counts_filtered)
}
