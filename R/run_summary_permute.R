source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)
library(tsutils)

forecast_horizon <- 1
source("R/metrics.R")

method_name <- function(representor, distance, cluster) {
  representor <- strsplit(representor, "-")[[1]][1]
  representor <- ifelse(!is.na(representor), switch(representor, error = "ER-", error.features = "ERF-",
                                                    ts = "TS-", ts.features = "TSF-"), "")
  distance <- ifelse(distance == "","",
                     switch(distance, euclidean="-EUC", dtw="-DTW")
  )
  cluster <- strsplit(cluster, "-")[[1]][1]
  cluster <- ifelse(
    !is.na(cluster),
    switch(cluster, natural="Natural", Kmedoids="ME", hcluster="HC", base="Base"),
    "Two-level"
  )
  paste0(representor, cluster, distance)
}


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
      polygon(rep(ranks.means[select], 4) + r.stat/2 * 
                c(1, 1, -1, -1), c(0, rep(cols.number + 1, 2), 
                                   0), col = "gray90", border = NA)
      points(ranks.means, 1:cols.number, pch = 20, lwd = 0.5)
      if (cols.number > 30) {
        target_loc <- which(names(ranks.means) == target)
        axis_at <- seq(target_loc, 1, -5)
        axis_at <- sort(axis_at[2:(length(axis_at))])
        axis_at <- c(axis_at, seq(target_loc, cols.number, 5))
        axis_labels <- paste0(labels[axis_at], 
                              " - ", 
                              sprintf("%1.2f", round(ranks.means[axis_at], 2)))
      }
      axis(2, at = (1:cols.number)[-axis_at], labels = FALSE)
      axis(2, at = axis_at, labels = axis_labels, las = 2, col.ticks="red", cex.axis=0.8)
      axis(1)
      for (i in 1:cols.number) {
        lines(ranks.means[i] + c(-1, 1) * 0.5 * r.stat, 
              rep(i, times = 2), type = "o", lwd = 1, col = pcol, 
              pch = 20)
      }
      idx <- abs(ranks.means[select] - ranks.means) < r.stat
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



mcb_series_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>% arrange(batch) %>%
    pull(rmsse) %>%
    do.call(c, .)
  rand_rmsse <- 
    rand %>% arrange(batch, cluster) %>%
    select(cluster, rmsse, batch) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    mutate_at("rmsse", purrr::map, function(g){
      g <- g %>% arrange(batch)
      list(do.call(c, g$rmsse))
    }) %>%
    tidyr::unnest(rmsse) %>%
    pull(rmsse) %>%
    do.call(cbind, .)
  all_rmsse <- cbind(orig_rmsse, rand_rmsse)
  colnames(all_rmsse) <- c(name, 1:100)
  nemenyi(all_rmsse, plottype = "vmcb", target = name)
}

mcb_hierarchy_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>% rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch) %>%
    pull(rmsse)
  
  rand_rmsse <- rand %>% rowwise() %>% 
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch, cluster) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    pull(rmsse) %>%
    lapply(function(x) { x %>% arrange(batch) %>% pull(rmsse) }) %>%
    do.call(cbind, .)
  
  all_rmsse <- cbind(orig_rmsse, rand_rmsse)
  colnames(all_rmsse) <- c(name, 1:100)
  nemenyi(all_rmsse, plottype = "vmcb", target = name)
}


rmsse_benchmarks <- NULL
rank_natural_tbl <- list()
rank_cluster_tbl <- list()

for (path in c("tourism", "mortality")) {
  dt <- readRDS(sprintf("%s/%s/eval.rds", path, "ets"))
  # MCB Test based on single series rmsse
  # pdf(sprintf("manuscript/figures/%s_natural_mcb_single.pdf", path),width = 8, height = 6)
  # natural_series <- mcb_series_rmsse(
  #   dt$dtb %>% filter(cluster == "natural"),
  #   dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
  #   "Natural"
  # )
  # dev.off()
  # which(names(natural_series$means) == "Natural")
  
  # MCB Test based on hierarchy rmsse
  pdf(sprintf("manuscript/figures/%s_natural_mcb_whole.pdf", path),width = 8, height = 6)
  natural_hierarchy <- mcb_hierarchy_rmsse(
    dt$dtb %>% filter(cluster == "natural"),
    dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
    "Natural"
  )
  dev.off()
  
  
  # base, natural, cluster rmsse
  bench_rmsse <- dt$dtb %>%
    filter(!startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(rmsse=mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster="base") %>%
            rowwise() %>% mutate(rmsse=mean(rmsse)))
  
  
  pdf(sprintf("manuscript/figures/%s_benchmark_mcb_whole.pdf", path),width = 8, height = 6, family = "Helvetica")
  par(mex=1.1)
  bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()
  
  
  clusters <- dt$dtb %>% select(cluster, distance, representor) %>%
    filter(!startsWith(cluster, "permute-")) %>%
    filter(distance != "") %>%
    unique()
  
  # cluster
  
  best_ <- bench_rmsse %>%
    group_by(representor, distance, cluster) %>%
    summarise(rmsse = mean(rmsse) * 100, .groups = "drop") %>%
    filter(representor != "") %>%
    arrange(rmsse)
  
  best_ <- best_[1,]
  
  best_name <- method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])
  
  # pdf(sprintf("manuscript/figures/%s_cluster_mcb_single.pdf", path),width = 8, height = 6)
  # cluster_series_rmsse <- mcb_series_rmsse(
  #   dt$dtb %>% filter(representor == best_$representor[[1]],
  #                     cluster == best_$cluster[[1]],
  #                     distance == best_$distance[[1]]),
  #   dt$dtb %>% filter(representor == best_$representor[[1]],
  #                     startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
  #                     distance == best_$distance[[1]]),
  #   best_name
  # )
  # dev.off()
  # cluster_series_rmsse$means
  
  
  pdf(sprintf("manuscript/figures/%s_cluster_mcb_whole.pdf", path),width = 8, height = 6)
  cluster_hierarchy_rmsse <- mcb_hierarchy_rmsse(
    dt$dtb %>% filter(representor == best_$representor[[1]],
                      cluster == best_$cluster[[1]],
                      distance == best_$distance[[1]]),
    dt$dtb %>% filter(representor == best_$representor[[1]],
                      startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
                      distance == best_$distance[[1]]),
    method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])
  )
  dev.off()
  
  bench_rmsse <- bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    group_by(method) %>%
    summarise(rmsse = mean(rmsse) * 100) %>%
    mutate(rmsse = round(rmsse, digits = 3)) %>%
    arrange(rmsse)
  
  if (path == "tourism") {
    rmsse_benchmarks <- bench_rmsse %>% rename(tourism = "rmsse")
  } else {
    rmsse_benchmarks <- bench_rmsse %>% 
      rename(mortality = "rmsse") %>%
      left_join(rmsse_benchmarks, by = "method")
  }
  rank_natural_tbl[[path]] <- c(
    bench_rmsse$rmsse[which(bench_rmsse$method == "Natural")],
    which(names(natural_hierarchy$means) == "Natural"))
  rank_cluster_tbl[[path]] <- c(
    bench_rmsse$rmsse[which(bench_rmsse$method == best_name)],
    which(names(cluster_hierarchy_rmsse$means) == best_name))
  
  if (path == "mortality") {
    # averaging
    rank_average_tbl <- list()
    all_rmsse3 <- NULL
    for (batch in unique(dt$dtb$batch)) {
      store_path <- sprintf("%s/%s/batch_%s.rds", path, "ets", batch)
      data <- readRDS(store_path)
      data_tibble <- nl2tibble(data$nl) %>% 
        filter(representor != "")
      
      data_tibble$permute <- sapply(data_tibble$cluster, function(x){
        if (!startsWith(x, "permute")) {
          return (0)
        }
        x <- strsplit(x, "-")[[1]]
        x <- x[length(x)]
        as.integer(x)
      })
      data_tibble$cluster <- sapply(data_tibble$cluster, function(x){
        if (!startsWith(x, "permute")) {
          return (x)
        }
        split_x <- strsplit(x, "-")[[1]]
        x <- stringi::stri_replace_all_fixed(x, "permute-", "")
        x <- stringi::stri_replace_all_fixed(x, paste0("-", split_x[length(split_x)]), "")
        x
      })
      
      avg <- data_tibble %>% select(representor, cluster, distance, rf, permute) %>%
        arrange(permute) %>%
        group_by(permute) %>%
        tidyr::nest(rf = -"permute") %>%
        mutate_at("rf", purrr::map, function(x) {
          sapply(x$rf, function(g) { g[["mint"]][1:forecast_horizon,,drop=FALSE] }, simplify = "array") %>%
            apply(c(1,2), mean)
        })
      if (is.null(dim(data$tts))){
        tts <- c(sum(data$tts), data$tts)
        tts <- matrix(tts, nrow=1)
      } else { 
        tts <- cbind(rowSums(data$tts), data$tts)[1:forecast_horizon,,drop=FALSE]
      }
      bts <- cbind(rowSums(data$bts), data$bts)
      
      rmsses <- sapply(iterators::iter(avg$rf), function(f) {
        sapply(1:NCOL(f), function(x){
          metric.rmsse(tts[,x], f[,x], bts[,x])
        }) %>% mean()
      })
      all_rmsse3 <- rbind(all_rmsse3, rmsses)
    }
    colnames(all_rmsse3) <- c("Average", 1:100)
    
    pdf("manuscript/figures/mortality_average_mcb_whole.pdf", height = 6, width = 8)
    test_ <- nemenyi(all_rmsse3, plottype = "vmcb", target = "Average")
    dev.off()
    rank_average_tbl[["mortality"]] <-
      c(round(colMeans(all_rmsse3)[1]*100, digits=3), which(names(test_$means) == "Average"))
  }
}


write.csv(rmsse_benchmarks, "manuscript/figures/cluster_rmsse.csv")

rank_natural_tbl <- data.frame(rank_natural_tbl)
row.names(rank_natural_tbl) <- c("RMSSE", "Rank")
write.csv(rank_natural_tbl, "manuscript/figures/rank_natural.csv")

rank_cluster_tbl <- data.frame(rank_cluster_tbl)
row.names(rank_cluster_tbl) <- c("RMSSE", "Rank")
write.csv(rank_cluster_tbl, "manuscript/figures/rank_cluster.csv")


rank_average_tbl <- data.frame(rank_average_tbl)
row.names(rank_average_tbl) <- c("RMSSE", "Rank")
write.csv(rank_average_tbl, "manuscript/figures/rank_average.csv")
