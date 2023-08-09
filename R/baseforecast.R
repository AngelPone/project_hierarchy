library(forecast)

f.arima <- function(x){
  mdl <- auto.arima(ts(x, frequency = 12))
  list(basef=forecast(mdl, h=12)$mean, resid=residuals(mdl))
}
