representator.ts <- function(x) {
  apply(x$bts, 2, function(x){
    (x - mean(x)) / sd(x)
  })
}

representator.accuracy <- function(x) {
  sapply(2:NCOL(x$resid), function(idx) {
    sqrt(mean(x$resid[,idx]^2) / mean(diff(x$bts[,idx-1], 12)^2))
  })
}

representator.forecast <- function(x) {
  apply(unclass(x$basef), 2, function(x){
    if (length(x) == 1) {
      return(0)
    }
    if (sd(x) < 1e-4) {
      return (vector("numeric", length(x)))
    }
    (x - mean(x)) / sd(x)
  })
}

representator.error <- function(x){
  stopifnot(!is.null(x$basef))
  n <- NROW(x$S)
  m <- NCOL(x$S)
  apply(x$resid[,(n-m+1):n], 2, function(x) {
    (x - mean(x)) / sd(x)
  })
}



features.compute <- function(data, frequency=frequency) {
  
  feature_lst <- c("acf_features", "arch_stat", "autocorr_features", "crossing_points", "dist_features",
                   "entropy", "heterogeneity", "hurst", "lumpiness", "stability", "pacf_features", "stl_features",
                   "unitroot_kpss", "unitroot_pp", "nonlinearity", "max_level_shift", "max_var_shift", "max_kl_shift",
                   "holt_parameters", "hw_parameters", "flat_spots")
  remove_zero_sd <- function(x){
    x[which(apply(x, 1, sd) > 0),]
  }
  
  data$features <- list()
  
  ts_features <- tsfeatures::tsfeatures(ts(data$bts, frequency = frequency), features = feature_lst)
  ts_features <- t(unname(as.matrix(ts_features[,!(colnames(ts_features) %in% c("nperiods", "seasonal_period"))])))
  ts_features <- remove_zero_sd(ts_features)
  
  error_features <- tsfeatures::tsfeatures(ts(data$resid[,2:NCOL(data$resid)], frequency = frequency), features = feature_lst)
  error_features <- t(unname(as.matrix(error_features[,!(colnames(error_features) %in% c("nperiods", "seasonal_period"))])))
  error_features <- remove_zero_sd(error_features)
  
  data$features$ts <- remove_zero_sd(ts_features)
  data$features$error <- error_features
  data
}


representator.ts.features <- function(x) {
  t(apply(x$features$ts, 1, function(g) {
    (g - mean(g)) / sd(g)
  }))
}

representator.error.features <- function(x) {
  t(apply(x$features$error, 1, function(g) {
    (g - mean(g)) / sd(g)
  }))
}
