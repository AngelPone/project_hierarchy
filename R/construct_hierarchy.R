library(forecast)
library(cluster)
library(foreach)
library(dplyr)
source("representator.R")
source("distance.R")
source("clustering.R")
source("reconciliation.R")
source("baseforecast.R")
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
      rf = NULL
    ),
    class = "hts"
  )
}



#' function to iterator over S, find the forecast store and return the forecast
#' 
fhts_helper <- function(S, all_ts, f_str, h, frequency) {
  f <- get(paste0("f.", f_str))
  stored_forecasts <- lapply(iterators::iter(S, by = "row"), function(row){
    f.store.read(row, f_str)
  })
  idx_to_forecast <- which(sapply(stored_forecasts, is.null))
  
  if (exists("num.cores")) {
    bf <- foreach::foreach(x = iterators::iter(all_ts[,idx_to_forecast,drop=FALSE], by = "column"), .packages = c("forecast")) %dopar% {
      f(x, h=h, frequency=frequency)
    }
  } else {
    bf <- apply(all_ts[,idx_to_forecast, drop=FALSE], 2, f, h=h, frequency=frequency)
  }
  for (i in seq_along(idx_to_forecast)) {
    S_idx <- idx_to_forecast[i]
    f.store.write(S[S_idx,], f_str, bf[[i]])
    stored_forecasts[[S_idx]] <- bf[[i]]
  }
  stored_forecasts
}

#' function to forecast hts
#' @param x hts
#' @param f_str baseforecast method string
forecast.hts <- function(x, f_str, h, frequency){
  if (is.null(x$basef)) {
    all_ts <- x$bts %*% t(x$S)
    
    bf <- fhts_helper(x$S, all_ts, f_str, h, frequency)
    
    x$basef <- unname(do.call(cbind, lapply(bf, function(x){x$basef})))
    x$resid <- unname(do.call(cbind, lapply(bf, function(x){x$resid})))
  }
  
  if (length(x$nl) > 0) {
    new_S <- do.call(rbind, lapply(x$nl, function(g){ g$S }))
    all_nts <- x$bts %*% t(new_S)
    bf <- fhts_helper(new_S, all_nts, f_str, h, frequency)
  }
  
  for (level in seq_along(x$nl)) {
    if (is.null(x$nl[[level]]$basef)) {
      all_nts <-x$bts %*% t(x$nl[[level]]$S)
      bf <- fhts_helper(x$nl[[level]]$S, all_nts, f_str, h, frequency)
      x$nl[[level]]$basef <- unname(do.call(cbind, lapply(bf, function(x){unclass(x$basef)})))
      
      if (is.null(dim(x$nl[[level]]$basef))) {
        x$nl[[level]]$basef <- matrix(x$nl[[level]]$basef, nrow = 1)
      }
      
      x$nl[[level]]$resid <- unname(do.call(cbind, lapply(bf, function(x){unclass(x$resid)})))
    }
  }
  x
}

#' add new level
updateLevel.hts <- function(x, new_S, meta) {
  append(x$nl, list(S = S, basef=NULL, rf=NULL, rf_nl=NULL, meta=meta))
}


is.hts <- function(x) {
  "hts" %in% class(x)
}


#' construct hierarchy
#' 
#' @param hts hts object
#' @param representator function to transform time series: 
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
    representator,
    distance,
    cluster,
    keep_old = FALSE,
    ...) {
  
  stopifnot(is.hts(hts))
  
  n <- NCOL(hts$bts)
  cluster_input <- representator(hts)
  distance_mat <- matrix(0, n, n)
  
  for(i in 1:n) {
    for (j in 1:i) {
      distance_mat[i, j] <- distance(cluster_input[, i], cluster_input[, j])
      distance_mat[j, i] <- distance_mat[i, j]
    }
  }
  
  group_result <- cluster(distance_mat, ...)
  
  # temporal function that group list to summing matrix
  nl <- lapply(group_result, function(x) {
    list(S = x, basef = NULL)
  })
  
  concat_str <- function(x){
    do.call(paste0, as.list(x))
  }
  
  # remove duplicated rows
  orig_rows <- apply(hts$S, 1, concat_str)
  if (keep_old) {
    orig_rows <- append(orig_rows, do.call(c, lapply(hts$nl, function(g) {apply(g$S, 1, concat_str)})))
  }
  
  for (i in seq_along(nl)) {
    nl[[i]]$S <-nl[[i]]$S[which(!(apply(nl[[i]]$S, 1, concat_str) %in% orig_rows)),,drop=FALSE]
  }
  
  nl <- Filter(function(x){NROW(x$S) > 0}, nl)
  
  
  if (keep_old) {
    hts$nl <- append(hts$nl, nl)
  } else {
    hts$nl <- nl
  }
  
  hts
}


#' #' lloyd k-means implementation
#' #' 
#' #' Standardization: each ts is divided by standard error of itself
#' #' Initialization: randomly choose n_clusters ts as initial centers
#' #' New centers: arithmetic mean of samples in cluster
#' #' Termination: distance between new center and old center smaller than 
#' #' tolerance, or exceed max iteration numbers
#' #' 
#' kmeans <- function(ts_mat, distance, n_clusters, max_iter = 100, tol=1e-4) {
#'   # standardization
#'   ts_mat <- apply(ts_mat, 2, function(x)( x/sd(x) ))
#'   
#'   n <- NCOL(ts_mat)
#'   # Initialization, forgy method
#'   center <- ts_mat[,sample(1:n, n_clusters)]
#'   
#'   iter_num <- 0
#'   
#'   while(TRUE) {
#'     iter_num <- iter_num + 1
#'     # compute new clusters
#'     grplst <- vector(mode = "list", length = n_clusters)
#'     
#'     minidx <- apply(ts_mat, 2, function(x){
#'       which.min(apply(center, 2, function(c) { distance(c, x) }))
#'     })
#'     
#'     for (i in seq_along(grplst)) {
#'       grplst[[i]] <- which(minidx == i)
#'     }
#'     
#'     # cluster new centers
#'     new_center <- 
#'       do.call(cbind, lapply(grplst, function(grp) { rowMeans(matrix(ts_mat[,grp], ncol = length(grp))) }))
#'     
#'     if ((distance(new_center, center) < tol) | (iter_num > max_iter)) {
#'       break
#'     } else {
#'       center <- new_center
#'     }
#'   }
#'   list(grplst)
#' }


