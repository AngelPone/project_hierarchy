library(forecast)

f.arima <- function(x){
  if(length(x) == 1) {
    stop("length == 1")
  }
  mdl <- auto.arima(ts(x, frequency = 12))
  list(basef=forecast(mdl, h=12)$mean, resid=residuals(mdl))
}

#' List to store the base forecast of aggregated time series
#' Base forecast of new aggregation will be saved into this store.


#' @param x: row of S matrix
#' @param method base forecast method
f.store.read <- function(x, method) {
  x_str <- do.call(paste0, as.list(x))
  if (!(x_str %in% names(BASEFORECAST_STORE[[method]]))) {
    return(NULL)
  }
  return(BASEFORECAST_STORE[[method]][[x_str]])
}

f.store.write <- function(x, method, fcasts) {
  x_str <- do.call(paste0, as.list(x))
  BASEFORECAST_STORE[[method]][[x_str]] <<- fcasts
}

