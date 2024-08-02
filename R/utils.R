library(forecast, quietly = TRUE)
library(cluster, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(furrr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(dtw, quietly = TRUE)

plan(multisession, workers=8)


source("expr_utils.R")

#' RMSSE
metric.rmsse <- function(obs, pred, hist) {
  sqrt(mean((obs - pred)^2) / mean(diff(hist, 12)^2))
}


#' kmedoids
#' 
cluster.kmedoids <- function(distance_mat, n_clusters) {
  
  # calculate silhouette and determine optimal number of clusters
  max.avgwidths <- 0
  max.n_cluster <- 1
  max.pr <- NULL
  for (n_cluster in n_clusters) {
    if (n_cluster == 1) next
    pr <- pam(distance_mat, k=n_cluster, diss=TRUE)
    clus_silwidth <- summary(silhouette(pr))
    # if (min(clus_silwidth$clus.avg.widths) < 0.1) next
    if (clus_silwidth$avg.width > max.avgwidths) {
      max.avgwidths <- clus_silwidth$avg.width
      max.n_cluster <- n_cluster
      max.pr <- pr
    }
  }
  if (max.n_cluster == 1) { return(NULL) }
  
  grp2S <- function(grp) {
    grpvec <- grp$clustering
    do.call(rbind, lapply(unique(grpvec), function(grp){
      S_row <- vector("numeric", NCOL(distance_mat))
      S_row[which(grpvec == grp)] <- 1
      S_row
    }))
  }
  
  list(S=grp2S(max.pr), 
       info = list(n_cluster = max.n_cluster, 
                   width = summary(silhouette(max.pr))$clus.avg.widths,
                   avgwidth = summary(silhouette(max.pr))$avg.width,
                   medoids = max.pr$id.med)
  )
}

#' hierarchical clustering
#' 
#' @param method linkage method: see ?cluster::agnes
cluster.hcluster <- function(distance_mat, method) {
  hc <- agnes(distance_mat, diss = FALSE, method = method,
              keep.data = FALSE, keep.diss = FALSE)
  
  S <- matrix(0, NROW(hc$merge) - 1, NROW(distance_mat))
  
  for (i in 1:NROW(S)) {
    cur_idx <- hc$merge[i,]
    for (idx in cur_idx) {
      if (idx < 0) {
        S[i, abs(idx)] <- 1
      } else {
        S[i, which(S[idx, ] == 1)] <- 1
      }
    }
  }
  list(S)
}


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

#' function to generate base forecasts and reconciled forecasts for all
#' new hierarchies
hts.nlf <- function(htst, h, frequency) {
  
  # only forecast hierarchies that did not produce reconciled forecasts
  idx2forecast <- which(sapply(htst$nl, function(x) {length(x$rf) == 0}))
  if (length(idx2forecast) == 0) {
    return(htst)
  }
  
  rfs <- future_map(htst$nl[idx2forecast], function(x) {
    # two-level hierarchy
    if (is.null(x$S)) {
      return(list(
        rf = reconcile.mint(rbind(rep(1, m), diag(m)), htst$basef, htst$resid)
      ))
    }
    mid_ts <- htst$bts %*% t(as.matrix(x$S))
    mid_forecasts <- 
      map(as.list(iterators::iter(mid_ts, by="column")), \(x) f.ets(x, h=h, frequency=frequency))
    S <- rbind(rep(1, m), as.matrix(x$S), diag(m))
    basef <- cbind(htst$basef[,1,drop=FALSE],
                   do.call(cbind, map(mid_forecasts, "basef")),
                   htst$basef[,2:NCOL(htst$basef)])
    resid <- cbind(htst$resid[,1,drop=FALSE],
                   do.call(cbind, map(mid_forecasts, "resid")),
                   htst$resid[,2:NCOL(htst$resid)])
    list(rf=reconcile.mint(S, basef, resid), basef=do.call(cbind, map(mid_forecasts, "basef")))
  })
  
  for (l in seq_along(idx2forecast)) {
    idxinnl <- idx2forecast[l]
    htst$nl[[idxinnl]]$rf <- rfs[[l]]$rf
    if (htst$nl[[idxinnl]]$cluster == "natural") {
      htst$nl[[idxinnl]]$other <- list(basef = rfs[[l]]$basef)
    }
  }
  htst
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




nl2tibble <- function(x) {
  tibble(
    representor = map_chr(x, "representor"),
    distance = map_chr(x, "distance"),
    cluster = map_chr(x, "cluster"),
    S = map(x, "S"),
    rf = map(x, "rf"),
    other = map(x, "other")
  )
}
