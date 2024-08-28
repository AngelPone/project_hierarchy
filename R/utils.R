library(forecast, quietly = TRUE)
library(cluster, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(furrr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(dtw, quietly = TRUE)
library(Matrix, quietly = TRUE)
plan(multisession, workers=8)


# some quantities basically used in all scripts
dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon + 1



#' RMSSE
metric.rmsse <- function(obs, pred, hist) {
  sqrt(mean((obs - pred)^2) / mean(diff(hist, 12)^2))
}


#' kmedoids clustering
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
  
  ts_features <- tsfeatures::tsfeatures(ts(data$bts, frequency = frequency), 
                                        features = feature_lst,
                                        parallel = TRUE)
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


#' function to produce mcb test plot
mcb_hierarchy_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch) %>%
    pull(rmsse)
  
  permute_idx <- sapply(strsplit(rand$cluster, "-"), function(x){
    as.integer(x[length(x)])
  })
  rand_rmsse <- rand %>% mutate(cluster=permute_idx) %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch, cluster) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    pull(rmsse) %>%
    lapply(\(x) x %>% arrange(batch) %>% pull(rmsse)) %>%
    do.call(cbind, .)
  
  all_rmsse <- cbind(orig_rmsse, rand_rmsse)
  colnames(all_rmsse) <- c(name, 1:100)
  nemenyi(all_rmsse, plottype = "vmcb", target = name)
}



method_name <- function(representor, distance, cluster) {
  if (startsWith(cluster, "permute")) {
    spl <- strsplit(cluster, "-")[[1]]
    return (spl[length(spl)])
  }
  representor <- strsplit(representor, "-")[[1]][1]
  representor <- ifelse(!is.na(representor), switch(representor,
                                                    error = "ER-",
                                                    error.features = "ERF-",
                                                    ts = "TS-",
                                                    ts.features = "TSF-"
  ), "")
  distance <- ifelse(distance == "", "",
                     switch(distance,
                            euclidean = "EUC",
                            dtw = "DTW"
                     )
  )
  cluster <- strsplit(cluster, "-")[[1]][1]
  cluster <- ifelse(
    !is.na(cluster),
    switch(cluster,
           natural = "Natural",
           Kmedoids = "-ME",
           hcluster = "-HC",
           base = "Base",
           "combination1" = "Combination",
           "combination2" = "Grouped"
    ),
    "Two-level"
  )
  paste0(representor, distance, cluster)
}


#' mcb test plot
#' original code from tsutils package
nemenyi <- function (data, conf.level = 0.95, sort = c(TRUE, FALSE), plottype = c("vline", 
                                                                                  "none", "mcb", "vmcb", "line", "matrix"), select = NULL, 
                     labels = NULL, target, ...) 
{
  sort <- sort[1]
  plottype <- match.arg(plottype, c("vline", "none", "mcb", 
                                    "vmcb", "line", "matrix"))
  if (length(dim(data)) != 2) {
    stop("Data must be organised as methods in columns and observations in rows.")
  }
  data <- as.matrix(data)
  data <- na.exclude(data)
  rows.number <- nrow(data)
  cols.number <- ncol(data)
  if (!is.null(select) && (select > cols.number)) {
    select <- NULL
  }
  if (plottype != "none") {
    sort <- TRUE
  }
  if (is.null(labels)) {
    labels <- colnames(data)
    if (is.null(labels)) {
      labels <- 1:cols.number
    }
  }
  else {
    labels <- labels[1:cols.number]
  }
  fried.pval <- stats::friedman.test(data)$p.value
  if (fried.pval <= 1 - conf.level) {
    fried.H <- "Ha: Different"
  }
  else {
    fried.H <- "H0: Identical"
  }
  r.stat <- stats::qtukey(conf.level, cols.number, Inf) * sqrt((cols.number * 
                                                                  (cols.number + 1))/(12 * rows.number))
  ranks.matrix <- t(apply(data, 1, function(x) {
    rank(x, na.last = "keep", ties.method = "average")
  }))
  ranks.means <- colMeans(ranks.matrix)
  ranks.intervals <- rbind(ranks.means - r.stat, ranks.means + 
                             r.stat)
  if (sort == TRUE) {
    order.idx <- order(ranks.means)
  }
  else {
    order.idx <- 1:cols.number
  }
  ranks.means <- ranks.means[order.idx]
  ranks.intervals <- ranks.intervals[, order.idx]
  labels <- labels[order.idx]
  if (!is.null(select)) {
    select <- which(order.idx == select)
  }
  if (plottype != "none") {
    args <- list(...)
    args.nms <- names(args)
    if (!("main" %in% args.nms)) {
      args$main <- paste0(
        sprintf("MCB test of %s dataset", args$dataset),
        "\nFriedman: ", format(round(fried.pval, 
                                     3), nsmall = 3), " (", fried.H, ") Critical distance: ", 
        format(round(r.stat, 3), nsmall = 3),sep = "")
      args$dataset <- NULL
    }
    if (!("xaxs" %in% names(args))) {
      args$xaxs <- "i"
    }
    if (!("yaxs" %in% names(args))) {
      args$yaxs <- "i"
    }
    nc <- max(nchar(labels))
    nc <- nc/1.75 + 1
    nr <- nchar(sprintf("%1.2f", round(max(ranks.means), 
                                       2)))/1.75
    parmar.def <- parmar <- graphics::par()$mar
  }
  if ((plottype == "mcb") | (plottype == "vmcb")) {
    cmp <- RColorBrewer::brewer.pal(3, "Set1")[1:2]
    if (fried.pval > 1 - conf.level) {
      pcol <- "gray"
    }
    else {
      pcol <- cmp[2]
    }
    mnmx <- range(ranks.means) + c(-0.5, 0.5) * r.stat
    mnmx <- mnmx + diff(mnmx) * 0.04 * c(-1, 1)
    if (plottype == "mcb") {
      if (!("xlab" %in% names(args))) {
        args$xlab <- ""
      }
      if (!("ylab" %in% names(args))) {
        args$ylab <- "Mean ranks"
      }
      if (is.null(args$xlim)) {
        args$xlim <- c(0, cols.number + 1)
      }
      if (is.null(args$ylim)) {
        args$ylim <- mnmx
      }
    }
    else {
      if (!("ylab" %in% names(args))) {
        args$ylab <- ""
      }
      if (!("xlab" %in% names(args))) {
        args$xlab <- "Mean ranks"
      }
      if (is.null(args$ylim)) {
        args$ylim <- c(0, cols.number + 1)
      }
      if (is.null(args$xlim)) {
        args$xlim <- mnmx
      }
    }
    args$x <- args$y <- NA
    args$axes <- FALSE
    if ((plottype == "mcb") && (parmar[1] < (nc + nr))) {
      parmar[1] <- nc + nr
    }
    if ((plottype == "vmcb") && (parmar[2] < (nc + nr))) {
      parmar[2] <- nc + nr
    }
    par(mar = parmar)
    if (is.null(select)) {
      select <- 1
    }
    do.call(plot, args)
    if (plottype == "mcb") {
      polygon(c(0, rep(cols.number + 1, 2), 0), rep(ranks.means[select], 
                                                    4) + r.stat/2 * c(1, 1, -1, -1), col = "gray90", 
              border = NA)
      points(1:cols.number, ranks.means, pch = 20, lwd = 10)
      axis(1, at = c(1:cols.number), labels = paste0(labels, 
                                                     " - ", sprintf("%1.2f", round(ranks.means, 2))), 
           las = 2)
      axis(2)
      for (i in 1:cols.number) {
        lines(rep(i, times = 2), ranks.means[i] + c(-1, 
                                                    1) * 0.5 * r.stat, type = "o", lwd = 1, col = pcol, 
              pch = 20)
      }
      idx <- abs(ranks.means[select] - ranks.means) < r.stat
      points((1:cols.number)[idx], ranks.means[idx], pch = 20, 
             lwd = 3, col = cmp[1])
    }
    else {
      target_loc <- which(names(ranks.means) == target)
      polygon(rep(ranks.means[target_loc], 4) + r.stat/2 * 
                c(1, 1, -1, -1), c(0, rep(cols.number + 1, 2), 
                                   0), col = "gray90", border = NA)
      points(ranks.means, 1:cols.number, pch = 20, lwd = 0.5)
      if (cols.number > 30) {
        axis_at <- seq(1, 101, 20)
        idx_at <- which(abs(target_loc - axis_at) <= 5)
        if (length(idx_at) == 1) {
          axis_at[idx_at] <- target_loc
        } else {
          axis_at <- c(axis_at, target_loc)
        }
        axis_labels <- paste0(labels[axis_at], 
                              " - ", 
                              sprintf("%1.2f", round(ranks.means[axis_at], 2)))
      }
      axis(2, at = (1:cols.number)[-axis_at], labels = FALSE)
      axis(2, at = axis_at, labels = axis_labels, las = 2, col.ticks="red", cex.axis=1.5)
      axis(1)
      for (i in 1:cols.number) {
        lines(ranks.means[i] + c(-1, 1) * 0.5 * r.stat, 
              rep(i, times = 2), type = "o", lwd = 1, col = pcol, 
              pch = 20)
      }
      idx <- abs(ranks.means[target_loc] - ranks.means) < r.stat
      points(ranks.means[idx], (1:cols.number)[idx], pch = 20, 
             lwd = 0.5, col = cmp[1])
    }
    box(which = "plot", col = "black")
  }
  if (plottype != "none") {
    par(mar = parmar.def)
  }
  return(structure(list(means = ranks.means, intervals = ranks.intervals, 
                        fpval = fried.pval, fH = fried.H, cd = r.stat, conf.level = conf.level, 
                        k = cols.number, n = rows.number), class = "nemenyi"))
}

