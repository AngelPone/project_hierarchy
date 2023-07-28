

#' hts
#' 
#' S: summing matrix
#' bts: bottom level series, should be a matrix
#' nl: new level
#" basef: base forecasts
#' rf: reconciled forecasts
hts <- function(S, bts) {
  structure(
    list(
      S = S, bts = bts,
      nl = list(),
      basef = NULL,
      rf = NULL
    ),
    class = "hts"
  )
}



#' function to forecast hts
forecast.hts <- function(x, f){
  if (is.null(x$basef)) {
    all_ts <- x$bts %*% t(x$S)
    x$basef <- apply(all_ts, 2, f)
  }
  if (!is.null(hts$nl$S) && is.null(hts$nl$basef)) {
    all_nts <-x$bts %*% t(x$nl$S)
    x$nl$basef <- apply(all_ts, 2, f)
  }
}

#' add new level
updateLevel.hts <- function(x, new_S, meta) {
  append(x$nl, list(S = S, basef=NULL, rf=NULL, rf_nl=NULL, meta=meta))
}

evaluate.hts <- function(x, metric = "RMSE") {
  if (is.null(x$rf)) {
    stop("Reconciled forecasts are null")
  }
  
  if (metric == "RMSE") {
    for (level in x$nl) {
      if (is.null(level$rf)) next
      level$acc <- sqrt(rowMeans((level$rf - x$bts)^2))
    }
  }
}


#' Euclidean distance
distance.euclidean <- function(x, y) { (x-y)^2 }


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
#' @param cluster function of clustering ( tts, distance ) -> group_lst
#' * tts is a matrix of transformed bts
#' * group_list should be a list of list. Outer list refers to multiple 
#' different clustering results, e.g, multiple run of K-means with different
#' cluter numbers, clustering path of hierarchical clustering.
#' The inner list refers to cluster member of each clusters.
#' The Outer list and inner list should both be indexed by integer.
#' @return hts object with new levels
build_level <- function(
    hts,
    representator,
    distance,
    cluster,
    keep_old = FALSE,
    ...) {
  
  stopifnot(is.hts(hts))
  
  cluster_input <- matrix(apply(hts$bts, 2, representator),
                         ncol = NCOL(hts$bts))
  group_result <- cluster(cluster_input, distance, ...)
  
  stopifnot(is.list(group_result[[1]]))
  
  # temporal function that group list to summing matrix
  grplst2Slst <- function(grplsts) {
    lapply(grplsts, function(grplst) {
      list(S = do.call(rbind, lapply(grp){
        S_row <- vector("numeric", NCOL(hts$bts))
        S_row[grp] <- 1
        S_row
      }), basef=NULL)
    })
  }
  
  nl <- grplst2Slst(group_result)
  
  if (keep_old) {
    hts$nl <- append(hts$nl, nl)
  } else {
    hts$nl <- nl
  }
  
  hts
}




#' lloyd k-means implementation
#' 
#' Standardization: each ts is divided by standard error of itself
#' Initialization: randomly choose n_clusters ts as initial centers
#' New centers: arithmetic mean of samples in cluster
#' Termination: distance between new center and old center smaller than 
#' tolerance, or exceed max iteration numbers
#' 
kmeans <- function(ts_mat, distance, n_clusters, max_iter = 100, tol=1e-4) {
  # standardization
  ts_mat <- apply(ts_mat, 2, function(x)( x/sd(x) ))
  
  n <- NCOL(ts_mat)
  # Initialization, forgy method
  center <- ts_mat[,sample(1:n, n_clusters)]
  
  iter_num <- 0
  
  while(TRUE) {
    iter_num <- iter_num + 1
    # compute new clusters
    grplst <- vector(mode = "list", length = n_clusters)
    
    minidx <- apply(ts_mat, 2, function(x){
      which.min(apply(center, 2, function(c) { distance(c, x) }))
    })
    
    for (i in seq_along(grplst)) {
      grplst[[i]] <- which(minidx == i)
    }
    
    # cluster new centers
    new_center <- 
      do.call(cbind, lapply(grplst, function(grp) { rowMeans(matrix(ts_mat[,grp], ncol = length(grp))) }))
    
    if ((distance(new_center, center) < tol) | (iter_num > max_iter)) {
      break
    } else {
      center <- new_center
    }
  }
  list(grplst)
}

# test_data <- NULL
# for (i in 1:100){
#   test_data <- rbind(test_data, rnorm(60, rep(c(10000, 20000, 30000, 40000, 50000, 100000), each=10)))
# }
# 
# debug(kmeans)
# kmeans(test_data, function(x, y){sqrt(sum((x-y)^2))}, 6)
