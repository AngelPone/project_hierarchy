library(forecast)

f.arima <- function(x){
  mdl <- auto.arima(x)
  list(basef=forecast(mdl, h=12)$mean, resid=residuals(mdl))
}
