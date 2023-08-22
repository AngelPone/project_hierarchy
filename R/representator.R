representator.ts <- function(x) {
  apply(x$bts, 2, function(x){
    x / sd(x)
  })
}


representator.error <- function(x){
  stopifnot(!is.null(x$basef))
  n <- NROW(x$S)
  m <- NCOL(x$S)
  apply(x$resid[,(n-m+1):n], 2, function(x) {
    x / sd(x)
  })
}

feature_lst <- c("acf_features", "arch_stat", "autocorr_features", "crossing_points", "dist_features",
                 "entropy", "heterogeneity", "hurst", "lumpiness", "stability", "pacf_features", "stl_features",
                 "unitroot_kpss", "unitroot_pp", "nonlinearity", "max_level_shift", "max_var_shift", "max_kl_shift",
                 "holt_parameters", "hw_parameters", "flat_spots")

features.compute <- function(data) {
  library(tsfeatures)
  FEATURES <<- list()
  ts_features <- tsfeatures(ts(data$bts, frequency = 12), features = feature_lst)
  FEATURES$ts <<- t(unname(as.matrix(ts_features[,!(colnames(ts_features) %in% c("nperiods", "seasonal_period"))])))
  
  error_features <- tsfeatures(ts(data$resid[,2:NCOL(data$resid)], frequency = 12), features = feature_lst)
  FEATURES$error <<- t(unname(as.matrix(error_features[,!(colnames(error_features) %in% c("nperiods", "seasonal_period"))])))
}


representator.ts.features <- function(x) {
  if (is.null(FEATURES)) {
    features.compute(x)
  }
  apply(FEATURES$ts, 1, function(g) {
    g / sd(g)
  })
}

representator.error.features <- function(x) {
  if (is.null(FEATURES)) {
    features.compute(x)
  }
  apply(FEATURES$error, 1, function(g) {
    g / sd(g)
  })
}