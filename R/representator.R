representator.ts <- function(x) {
  apply(x$bts, 2, function(x){
    (x - mean(x)) / sd(x)
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

feature_lst <- c("acf_features", "arch_stat", "autocorr_features", "crossing_points", "dist_features",
                 "entropy", "heterogeneity", "hurst", "lumpiness", "stability", "pacf_features", "stl_features",
                 "unitroot_kpss", "unitroot_pp", "nonlinearity", "max_level_shift", "max_var_shift", "max_kl_shift",
                 "holt_parameters", "hw_parameters", "flat_spots")

features.compute <- function(data, frequency=frequency) {
  library(tsfeatures)
  FEATURES <<- list()
  ts_features <- tsfeatures(ts(data$bts, frequency = frequency), features = feature_lst)
  FEATURES$ts <<- t(unname(as.matrix(ts_features[,!(colnames(ts_features) %in% c("nperiods", "seasonal_period"))])))
  
  remove_zero_sd <- function(x){
    x[which(apply(x, 1, sd) > 0),]
  }
  
  FEATURES$ts <<- remove_zero_sd(FEATURES$ts)
  
  error_features <- tsfeatures(ts(data$resid[,2:NCOL(data$resid)], frequency = frequency), features = feature_lst)
  FEATURES$error <<- t(unname(as.matrix(error_features[,!(colnames(error_features) %in% c("nperiods", "seasonal_period"))])))
  FEATURES$error <<- remove_zero_sd(FEATURES$error)
}


representator.ts.features <- function(x) {
  if (is.null(FEATURES)) {
    stop("FEATURES should be computed first!")
  }
  t(apply(FEATURES$ts, 1, function(g) {
    (g - mean(g)) / sd(g)
  }))
}

representator.error.features <- function(x) {
  if (is.null(FEATURES)) {
    stop("FEATURES should be computed first!")
  }
  t(apply(FEATURES$error, 1, function(g) {
    (g - mean(g)) / sd(g)
  }))
}
