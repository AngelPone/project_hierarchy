args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]


dt <- readRDS(sprintf("%s/%s/eval.rds", path, bfmethod))
dt_orig <- readRDS(sprintf("%s/%s/batch_0.rds", path, bfmethod))
dt_names <- readRDS(sprintf("%s/names.rds", path))
source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)



nemenyi <- function (data, conf.level = 0.95, sort = c(TRUE, FALSE), plottype = c("vline", 
                                                                       "none", "mcb", "vmcb", "line", "matrix"), select = NULL, 
          labels = NULL, ...) 
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
      points(1:cols.number, ranks.means, pch = 20, lwd = 3)
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
      points(ranks.means, 1:cols.number, pch = 20, lwd = 3)
      axis(2, at = c(1:cols.number), labels = paste0(labels, 
                                                     " - ", sprintf("%1.2f", round(ranks.means, 2))), 
           las = 2)
      axis(1)
      for (i in 1:cols.number) {
        lines(ranks.means[i] + c(-1, 1) * 0.5 * r.stat, 
              rep(i, times = 2), type = "o", lwd = 1, col = pcol, 
              pch = 20)
      }
      idx <- abs(ranks.means[select] - ranks.means) < r.stat
      points(ranks.means[idx], (1:cols.number)[idx], pch = 20, 
             lwd = 3, col = cmp[1])
    }
    box(which = "plot", col = "black")
  }
  if (plottype == "matrix") {
    rline <- array(NA, c(cols.number, 2))
    nem.mat <- array(0, c(cols.number, cols.number))
    for (i in 1:cols.number) {
      rline[i, ] <- c(which(ranks.means > ranks.intervals[1, 
                                                          i])[1], tail(which(ranks.means < ranks.intervals[2, 
                                                                                                           i]), 1))
      nem.mat[i, rline[i, 1]:rline[i, 2]] <- 1
    }
    diag(nem.mat) <- 2
    if (parmar[1] < nc) {
      parmar[1] <- nc
    }
    nr <- nchar(sprintf("%1.2f", round(max(ranks.means), 
                                       2)))/1.75
    if (parmar[2] < (nc + nr)) {
      parmar[2] <- nc + nr
    }
    par(mar = parmar)
    if (!("xlab" %in% names(args))) {
      args$xlab <- ""
    }
    if (!("ylab" %in% names(args))) {
      args$ylab <- ""
    }
    args$x <- 1:cols.number
    args$y <- 1:cols.number
    args$z <- nem.mat[order(order.idx), cols.number:1]
    args$axes <- FALSE
    cmp <- c("white", RColorBrewer::brewer.pal(3, "Set1")[2], 
             "black")
    if (fried.pval > 1 - conf.level) {
      cmp[3] <- "gray60"
      cmp[2] <- "gray70"
    }
    args$col <- cmp
    do.call(image, args)
    for (i in 1:cols.number) {
      abline(v = i + 0.5)
      abline(h = i + 0.5)
    }
    axis(1, at = 1:cols.number, labels = labels[order(order.idx)], 
         las = 2)
    axis(2, at = 1:cols.number, labels = paste0(rev(labels), 
                                                " - ", sprintf("%1.2f", round(rev(ranks.means), 2))), 
         las = 2)
    box()
    if (!is.null(select)) {
      polygon(order.idx[select] + c(-0.5, 0.5, 0.5, -0.5), 
              which((nem.mat[order(order.idx), cols.number:1])[order.idx[select], 
              ] == 2) + c(-0.5, -0.5, 0.5, 0.5), col = "red")
    }
  }
  if ((plottype == "line") | (plottype == "vline")) {
    rline <- matrix(NA, nrow = cols.number, ncol = 2)
    for (i in 1:cols.number) {
      tloc <- which((abs(ranks.means - ranks.means[i]) < 
                       r.stat) == TRUE)
      rline[i, ] <- c(min(tloc), max(tloc))
    }
    rline <- unique(rline)
    rline <- rline[apply(rline, 1, min) != apply(rline, 1, 
                                                 max), ]
    if (length(rline) == 2) {
      rline <- as.matrix(rline)
      rline <- t(rline)
    }
    k <- nrow(rline)
    cmp <- colorRampPalette(RColorBrewer::brewer.pal(12, 
                                                     "Paired"))(k)
    if (fried.pval > 1 - conf.level) {
      cmp <- rep("gray", times = k)
    }
    lbl <- paste0(labels, " - ", sprintf("%1.2f", round(ranks.means, 
                                                        2)))
    if (!("ylab" %in% names(args))) {
      args$ylab <- ""
    }
    if (!("xlab" %in% names(args))) {
      args$xlab <- ""
    }
    args$x <- args$y <- NA
    args$axes <- FALSE
    if (plottype == "line") {
      if (is.null(args$xlim)) {
        args$xlim <- c(1, cols.number)
      }
      if (is.null(args$ylim)) {
        args$ylim <- c(0, k + 1)
      }
    }
    else {
      if (is.null(args$xlim)) {
        args$xlim <- c(0, k + 1)
      }
      if (is.null(args$ylim)) {
        args$ylim <- c(1, cols.number)
      }
    }
    if ((plottype == "line") && (parmar[1] < (nc + nr))) {
      parmar[1] <- nc + nr
    }
    if ((plottype == "vline") && (parmar[2] < (nc + nr))) {
      parmar[2] <- nc + nr
    }
    par(mar = parmar)
    do.call(plot, args)
    if (plottype == "line") {
      points(1:cols.number, rep(0, cols.number), pch = 20, 
             lwd = 4)
      if (k > 0) {
        for (i in 1:k) {
          lines(rline[i, ], c(i, i), col = cmp[i], lwd = 4)
          lines(rep(rline[i, 1], times = 2), c(0, i), 
                col = "gray", lty = 2)
          lines(rep(rline[i, 2], times = 2), c(0, i), 
                col = "gray", lty = 2)
        }
      }
      axis(1, at = c(1:cols.number), labels = lbl, las = 2)
      if (!is.null(select)) {
        points(select, 0, pch = 20, col = RColorBrewer::brewer.pal(3, 
                                                                   "Set1")[1], cex = 2)
      }
    }
    else {
      points(rep(0, cols.number), 1:cols.number, pch = 20, 
             lwd = 4)
      if (k > 0) {
        for (i in 1:k) {
          lines(c(i, i), rline[i, ], col = cmp[i], lwd = 4)
          lines(c(0, i), rep(rline[i, 1], times = 2), 
                col = "gray", lty = 2)
          lines(c(0, i), rep(rline[i, 2], times = 2), 
                col = "gray", lty = 2)
        }
      }
      axis(2, at = c(1:cols.number), labels = lbl, las = 2)
      if (!is.null(select)) {
        points(0, select, pch = 20, col = RColorBrewer::brewer.pal(3, 
                                                                   "Set1")[1], cex = 2)
      }
    }
  }
  if (plottype != "none") {
    par(mar = parmar.def)
  }
  return(structure(list(means = ranks.means, intervals = ranks.intervals, 
                        fpval = fried.pval, fH = fried.H, cd = r.stat, conf.level = conf.level, 
                        k = cols.number, n = rows.number), class = "nemenyi"))
}

rank_compare_summary <- function(dt) {
  
  methods_random <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(cluster %in% c("random-15-10", "random-natural-10")) %>% 
    unique()
  methods_random$name <- ifelse(methods_random$cluster == "random-15-10", "FC-R", "FC-N")
  
  methods_dr <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(endsWith(cluster, "-dr")) %>% unique()
  
  switch_name <- function(x){
    x <- strsplit(x, "-")[[1]][1]
    switch(x, ts="TS", ts.features="TSF", error="ER", error.features = "ERF")
  }
  
  methods_dr$name <- paste0(sapply(methods_dr$representor, switch_name), "-",
                            ifelse(methods_dr$cluster == "hcluster-dr", "HC", "ME"))
  methods_dtw <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(distance == "dtw", representor %in% c("ts", "error")) %>% unique()
  methods_dtw$name <- paste0(ifelse(methods_dtw$representor == "ts", "TS", "ER"), "-",
                             ifelse(methods_dtw$cluster == "hcluster", "HC", "ME"), "-DTW")
  methods_natural <- dt$dtb %>% select(representor, distance, cluster) %>%
    filter(cluster %in% c("base", "", "natural", "cluster-average")) %>% unique()
  methods_natural$name <- ifelse(methods_natural$cluster == "base", "BASE",
                                 ifelse(methods_natural$cluster == "", "Original",
                                        ifelse(methods_natural$cluster == "natural", "Natural", "FC-C")))
  
  methods_df <- rbind(methods_random, methods_dr, methods_dtw, methods_natural)

  dtb <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    select(representor, distance, cluster, rmse, batch, name) %>%
    tidyr::nest(values = c("batch", "rmse")) %>%
    mutate_at("values", purrr::map, function(g) {
      g <- arrange(g, "batch") %>% pull("rmse")
      do.call(c, lapply(g, function(x) x[["mint"]]))
    })

  metric_mat <- do.call(cbind, dtb$values)
  dt_base <- do.call(c, dt$base[["rmse"]])
  metric_mat <- cbind(metric_mat, dt_base)
  
  colnames(metric_mat) <- c(dtb$name, "Base")
  nemenyi(metric_mat, plot = "vmcb", dataset = path)
}


pdf(sprintf("manuscript/figures/%s_mcb.pdf", path), width = 10, height = 6)
par(mar = c(5.1, 10, 4.1, 2.1))
rank_compare_summary(dt)
dev.off()





rank_compare_number <- function(dt) {
  
  methods_random <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(cluster %in% paste0(rep(c("random-15-", "random-natural-"), 3), rep(c(10, 20, 50), each=2))) %>% 
    unique()
  name <- ifelse(startsWith(methods_random$cluster, "random-15"), "FC-R", "FC-N")
  name2 <- sapply(methods_random$cluster, function(x){
    paste0("-", strsplit(x, "-")[[1]][3])
  })
  name2[name2=="-10"] <- ""
  methods_random$name <- paste0(name, name2)

  methods_natural <- dt$dtb %>% select(representor, distance, cluster) %>%
    filter(cluster == "cluster-average") %>% unique()
  methods_natural$name <- "FC-C"
  
  methods_df <- rbind(methods_random, methods_natural)
  
  dtb <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    select(representor, distance, cluster, rmse, batch, name) %>%
    tidyr::nest(values = c("batch", "rmse")) %>%
    mutate_at("values", purrr::map, function(g) {
      g <- arrange(g, "batch") %>% pull("rmse")
      do.call(c, lapply(g, function(x) x[["mint"]]))
    })
  
  metric_mat <- do.call(cbind, dtb$values)
  colnames(metric_mat) <- dtb$name
  
  pdf(sprintf("manuscript/figures/%s_number_mcb.pdf", path), width = 6, height = 6)
  par(mar = c(5.1, 10, 4.1, 2.1))
  a <- tsutils::nemenyi(metric_mat, plot = "vmcb")
  dev.off()
  
  # number of nodes
  n_repeats <- sapply(strsplit(methods_df$name, "-"), function(x) {
    if (x[[2]] == 'C') return(1)
    ifelse(length(x) == 3, as.integer(x[[3]]), 10) 
  })

  natural_node <- dt$dtb %>% filter(cluster=='natural') %>% pull(S) %>% sapply(NROW) %>% mean()
  method <- sapply(strsplit(methods_df$name, "-"), function(x) {x[[2]]})
  cluster_node <- dt$dtb %>% rowwise() %>% 
    mutate(n_node = NROW(S)) %>% 
    group_by(batch) %>% summarise(n_node = sum(n_node)) %>% 
    pull(n_node) %>% mean()
  method_node <- ifelse(method == "N", natural_node, ifelse(method == "R", 15, as.integer(cluster_node)))
  methods_df <- methods_df %>% mutate(n_node = method_node * n_repeats)
  methods_df %>% mutate(ranks = a$means[methods_df$name]) %>% arrange(desc(ranks))
}


df <- rank_compare_number(dt)
pdf(sprintf("manuscript/figures/%s_number_node.pdf", path), width = 6, height = 5.5)
par(mar = c(5.1, 10, 4.1, 2.1))

ggplot(data = df, mapping = aes(x = factor(name, levels = rev(name)), y=n_node)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  xlab("") +
  ylab("Number of new series") +
  theme_grey(base_size = 15)
dev.off()




