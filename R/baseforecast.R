library(forecast, quietly = TRUE)

f.arima <- function(x, h, frequency){
  mdl <- auto.arima(ts(x, frequency = frequency))
  list(basef=forecast(mdl, h=h)$mean, resid=residuals(mdl, type = "response"))
}

f.ets <- function(x, h, frequency) {
  mdl <- ets(ts(x, frequency = frequency))
  list(basef=forecast(mdl, h=h)$mean, resid=residuals(mdl, type = "response"))
}
