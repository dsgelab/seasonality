source('seasonality_dashboard/utils.R')
a <- 13;b <- 1
k_trend <- 6
k_seasonal <- 6
dat$nr_days <- 30
f_null <- 'COUNT ~ s(EVENT_YEAR, k=k_trend, bs="ps")'
offset <- 'offset(log(nr_days))'
f_null <- paste(f_null,offset,sep='+')
f_seasonal <- paste(f_null,'s(EVENT_MONTH, k=k_seasonal, bs="cp")',sep='+')
mod_brms <- brm(bf(as.formula(f_seasonal)),
                      data = dat, family = negbinomial(), cores = 4, seed = 17,
                      knots=list(EVENT_MONTH = c(a-0.5, b-0.5)),
                      iter = 4000, warmup = 1000, thin = 10, refresh = 0,
                      control = list(adapt_delta = 0.99))

mod_adj <- run_seasonality_gam(dat,a=1,b=13,seasonal_spline_avg = 'hax')

plot(conditional_smooths(mod_brms,smooths='s(EVENT_MONTH, k=k_seasonal, bs="cp")'))
plot(conditional_effects(mod_brms,effects='EVENT_MONTH'))
