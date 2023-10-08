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

nl2tibble <- function(x) {
  output <- vector("list", 6)
  names(output) <- c("representor", "cluster", "distance", "S", "rf", "other")
  
  for (i in seq_along(x)) {
    for (n in names(x[[i]])) {
      if (is.character(x[[i]][[n]])) {
        output[[n]] <- append(output[[n]], x[[i]][[n]])
      } else {
        output[[n]] <- append(output[[n]], list(x[[i]][[n]]))
      }
    }
  }
  tibble::as_tibble(output)
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
    representor,
    distance,
    cluster,
    ...) {
  
  stopifnot(is.hts(hts))
  
  n <- NCOL(hts$bts)
  
  distance_mat <- hts$distance[[representor]][[distance]]
  
  cluster(distance_mat, ...)
}



