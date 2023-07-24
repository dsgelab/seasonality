get_pred_gamm <- function(mod,dat,type='brms',tr=function(x) x){
  if(type=='brms'){
    pred_dat <- predict(mod,newdata=dat) %>%
      as_tibble() %>%
      select(mean=Estimate,lower=Q2.5,upper=Q97.5)
  }else{
    pred_est <- tr(predict(mod,new_data=dat))
    pred_dat <- tibble(mean=pred_est,lower=pred_est,upper=pred_est)
  }
  return(pred_dat)
}

theme_font <- theme(plot.title = element_text(size=18),
                    strip.text = element_text(size=16),
                    axis.text.x = element_text(size=16),
                    axis.text.y = element_text(size=16),
                    axis.title.x = element_text(size=16),
                    axis.title.y = element_text(size=16),
                    legend.title = element_text(size=18),
                    legend.text = element_text(size=16))
