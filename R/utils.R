library(forecast, quietly = TRUE)
library(cluster, quietly = TRUE)
library(foreach, quietly = TRUE)
library(dplyr, quietly = TRUE)

source("clustering.R")
source("metrics.R")
source("expr_utils.R")


#' hts
#' 
#' S: summing matrix
#' bts: bottom level series, should be a matrix
#' nl: new level
#" basef: base forecasts
#' rf: reconciled forecasts
hts <- function(S, bts, tts) {
  structure(
    list(
      S = S, bts = bts, tts = tts,
      nl = list(),
      basef = NULL,
      resid = NULL,
      features = NULL,
      distance = NULL
    ),
    class = "hts"
  )
}


add_nl <- function(data, S, representor, distance, cluster, other=NULL) {
  if (!is.null(S)) { S <- as(S, "sparseMatrix") }
  data$nl[[length(data$nl)+1]] <- list(representor = representor, distance = distance, 
                       cluster = cluster, other = other,
                       S=S, rf=list())
  data
}

library(forecast)
library(furrr)
library(purrr)
#' function to produce base forecast for single time series using ets model
f.ets <- function(x, h, frequency) {
  mdl <- ets(ts(x, frequency = frequency))
  list(basef=as.numeric(forecast(mdl, h=h)$mean), 
       resid=as.numeric(residuals(mdl, type = "response")))
}

#' function to produce base forecast for hts using ets model
#' @param x hts
#' @param h forecast horizon
#' @param frequency frequency
hts.basef <- function(x, h, frequency) {
  all_ts <- x$bts %*% t(x$S)
  bf <- future_map(as.list(iterators::iter(all_ts, by = "column")),
                   \(x) f.ets(x, h=h, frequency))
  x$basef <- unname(do.call(cbind, map(bf, "basef")))
  x$resid <- unname(do.call(cbind, map(bf, "resid")))
  x
}

hts.nlf <- function(htst, f_str, h, frequency) {
  f <- get(paste0("f.", f_str))
  
  idx2forecast <- which(sapply(htst$nl, function(x) {length(x$rf) == 0}))
  if (length(idx2forecast) == 0) {
    return(htst)
  }
  smat2sstr <- function(S) {
    S_str <- c()
    for (i in 1:NROW(S)) {
      S_str <- c(S_str, do.call(paste0, as.list(S[i,])))
    }
    S_str
  }
  
  S <- do.call(rbind, lapply(htst$nl[idx2forecast], function(g) { g$S }))
  print(sprintf("totally %s series are constructed", NROW(S)))
  S_str <- smat2sstr(as.matrix(S))
  fcasts <- list()
  
  for (idx in seq_along(S_str)) {
    if (S_str[idx] %in% names(fcasts)) {
      next
    } else {
      fcasts[[S_str[idx]]] <- S[idx, ]
    }
  }
  
  S_toforecast <- do.call(rbind, fcasts)
  print(sprintf("totally %s series are forecast", NROW(S_toforecast)))
  
  allts <- htst$bts %*% t(S_toforecast)
  
  bf <- foreach::foreach(x = iterators::iter(allts, by = "column"), .packages = c("forecast")) %dopar% {
    f(x, h=h, frequency=frequency)
  }
  names(bf) <- names(fcasts)
  
  rfs <- foreach(nl=iterators::iter(htst$nl[idx2forecast])) %do% {
    
    S_nl <- NULL
    if (!is.null(nl$S)) {
      S_nl <- as.matrix(nl$S)
    }
    
    S_str <- smat2sstr(S_nl)
    basef_nl <- do.call(cbind, lapply(S_str, function(g) { bf[[g]]$basef }))
    resid_nl <- do.call(cbind, lapply(S_str, function(g) { bf[[g]]$resid }))
    
    basef <- cbind(htst$basef[,1,drop=FALSE], 
                   basef_nl, 
                   htst$basef[,2:NCOL(htst$basef),drop=FALSE])
    resid <- cbind(htst$resid[,1], resid_nl, htst$resid[,2:NCOL(htst$basef)])
    

    S <- rbind(rep(1, NCOL(data$S)), S_nl, diag(NCOL(data$S)))
    
    reconcile.all(S, basef, resid)
  }
  
  for (l in seq_along(idx2forecast)) {
    idxinnl <- idx2forecast[l]
    htst$nl[[idxinnl]]$rf <- rfs[[l]]
  }
  htst
}


is.hts <- function(x) {
  "hts" %in% class(x)
}


#' construct hierarchy
#' 
#' @param hts hts object
#' @param representor function to transform time series: 
#' ( ts ) -> transformed_ts
#' * With same length time series input, representation should return same
#' length time series/point output.
#' @param distance function of similarity/diversity measure:
#' ( tts1, tts2 ) -> Distance(tts1, tts2)
#' @param cluster function of clustering ( distance, ... ) -> group_lst
#' * distance is the distance matrix
#' * group_list should be a list. Outer list refers to multiple 
#' different clustering results, e.g, multiple run of K-means with different
#' cluter numbers, clustering path of hierarchical clustering.
#' Each list is a vector indicating which cluster each series belongs to
#' @return hts object with new levels
build_level <- function(
    hts,
    representor,
    distance,
    cluster,
    ...) {
  
  stopifnot(is.hts(hts))
  
  n <- NCOL(hts$bts)
  
  distance_mat <- hts$distance[[representor]][[distance]]
  
  cluster(distance_mat, ...)
}



#' time series representation 
representor.ts <- function(x, dr = FALSE) {
  res <- apply(x$bts, 2, function(x){
    (x - mean(x)) / sd(x)
  })
  if (dr) {
    res <- prcomp(t(res), scale.=TRUE)
    vari <- which(cumsum((res$sdev^2)/sum(res$sdev^2)) > 0.8)
    vari <- min(vari, 10)
    res <- t(res$x[,1:vari])
  }
  res
}

#' time series representations with error
representor.error <- function(x, dr = FALSE){
  stopifnot(!is.null(x$basef))
  n <- NROW(x$S)
  m <- NCOL(x$S)
  res <- apply(x$resid[,(n-m+1):n], 2, function(x) {
    (x - mean(x)) / sd(x)
  })
  if (dr) {
    res <- prcomp(t(res), scale.=TRUE)
    vari <- which(cumsum((res$sdev^2)/sum(res$sdev^2)) > 0.8)
    vari <- min(vari, 10)
    res <- t(res$x[,1:vari])
  }
  res
}


#' function to compute features
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


representor.ts.features <- function(x, dr=FALSE) {
  res <- t(apply(x$features$ts, 1, function(g) {
    (g - mean(g)) / sd(g)
  }))
  if (dr) {
    res <- prcomp(t(res), scale.=TRUE)
    vari <- which(cumsum((res$sdev^2)/sum(res$sdev^2)) > 0.8)
    vari <- min(vari, 10)
    res <- t(res$x[,1:vari])
  }
  res
}

representor.error.features <- function(x, dr=FALSE) {
  res <- t(apply(x$features$error, 1, function(g) {
    (g - mean(g)) / sd(g)
  }))
  if (dr) {
    res <- prcomp(t(res), scale.=TRUE)
    vari <- which(cumsum((res$sdev^2)/sum(res$sdev^2)) > 0.8)
    vari <- min(vari, 10)
    res <- t(res$x[,1:vari])
  }
  res
}

#' distance
library(dtw, quietly = TRUE)
#' euclidean distance
distance.euclidean <- function(x, y) { sqrt(sum((x - y)^2)) }
#' dtw distance
distance.dtw <- function(x, y) { dtw(x, y, distance.only = TRUE)$distance }





# Reconciliation functions
reconcile.mint <- function(S, basef, resid){
  idx <- c(which(rowSums(S) > 1), (NROW(S) - NCOL(S) + 1):NROW(S))
  C <- S[which(rowSums(S) > 1),,drop=FALSE]
  basef <- basef[,idx]
  resid <- resid[,idx]
  n <- length(idx)
  m <- NCOL(S)
  unname(FoReco::csrec(basef, comb="shr", agg_mat=C, res=resid)[, c(1, (n-m+1):n)])
}

