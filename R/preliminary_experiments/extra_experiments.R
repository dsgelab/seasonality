
endpoint_nr <- 1
monthly_counts <-filter(monthly_counts_FinRegistry,ENDPOINT_NR==endpoint_nr)
monthly_counts$pred <- exp(predict(holiday_adj_gam$gam,newdata=expand_grid(EVENT_YEAR=seq(1998,2019),EVENT_MONTH=seq_len(12)) %>% mutate(nr_days=30,ENDPOINT_NR=endpoint_nr)))
ggplot(mutate(monthly_counts,month=1:n())) +
  geom_point(aes(x=EVENT_DATE,
                 y=COUNT)) +
  geom_line(aes(x=EVENT_DATE,
                y=pred))

# monthly_counts_FinRegistry$nr_days <- 30
# monthly_counts_FinRegistry <- mutate(monthly_counts_FinRegistry,nr_days=30) %>%
#                               group_by(ENDPOINT) %>%
#                               mutate(ENDPOINT_NR=cur_group_id()) %>%
#                               ungroup() %>%
#                               filter(ENDPOINT_NR<10)
# holiday_adj_gam <- gamm(COUNT ~ offset(log(nr_days)) +
#                                s(EVENT_YEAR,ENDPOINT_NR,bs='fs',k=11,xt=list(bs="tp")) +
#                                s(EVENT_MONTH,ENDPOINT_NR,bs='fs',k=6,xt=list(bs="tp")),
#                         random=list(ENDPOINT_NR=~1),
#                        family=quasipoisson(), data=monthly_counts_FinRegistry)
