library(forecast, quietly = TRUE)
library(cluster, quietly = TRUE)
library(foreach, quietly = TRUE)
library(dplyr, quietly = TRUE)
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

hts.basef <- function(x, f_str, h, frequency) {
  if (is.null(x$basef)) {
    f <- get(paste0("f.", f_str))
    all_ts <- x$bts %*% t(x$S)
    bf <- foreach::foreach(x = iterators::iter(all_ts, by = "column"), .packages = c("forecast")) %dopar% {
      f(x, h=h, frequency=frequency)
    }
    x$basef <- unname(do.call(cbind, lapply(bf, function(x){x$basef})))
    x$resid <- unname(do.call(cbind, lapply(bf, function(x){x$resid})))
  }
  x
}

hts.nlf <- function(htst, f_str, h, frequency) {
  f <- get(paste0("f.", f_str))
  
  
  smat2sstr <- function(S) {
    S_str <- c()
    for (i in 1:NROW(S)) {
      S_str <- c(S_str, do.call(paste0, as.list(S[i,])))
    }
    S_str
  }
  
  S <- do.call(rbind, lapply(htst$nl, function(g) { g$S }))
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
  
  rfs <- foreach(nl=iterators::iter(htst$nl)) %do% {
    
    S_nl <- NULL
    if (!is.null(nl$S)) {
      S_nl <- as.matrix(nl$S)
    }
    
    S_str <- smat2sstr(S_nl)
    basef_nl <- do.call(cbind, lapply(S_str, function(g) { bf[[g]]$basef }))
    resid_nl <- do.call(cbind, lapply(S_str, function(g) { bf[[g]]$resid }))
    
    basef <- cbind(htst$basef[,1], basef_nl, htst$basef[,2:NCOL(htst$basef)])
    resid <- cbind(htst$resid[,1], resid_nl, htst$resid[,2:NCOL(htst$basef)])
    

    S <- rbind(rep(1, NCOL(data$S)), S_nl, diag(NCOL(data$S)))
    
    reconcile.all(S, basef, resid)
  }
  
  for (l in seq_along(htst$nl)) {
    htst$nl[[l]]$rf <- rfs[[l]]
  }
  htst
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
    representor,
    distance,
    cluster,
    ...) {
  
  stopifnot(is.hts(hts))
  
  n <- NCOL(hts$bts)
  
  distance_mat <- hts$distance[[representor]][[distance]]
  
  cluster(distance_mat, ...)
}



