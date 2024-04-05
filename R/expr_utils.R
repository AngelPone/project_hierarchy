library(Matrix, quietly = TRUE)

create_new_output <- function() {
  output <- vector("list", 5)
  names(output) <- c("representator", "distance", "cluster", "accuracy", "other")
  output
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


add_result <- function(output, representotar, distance, 
                       cluster, rf, other = NULL) {
  if (is.null(output$representator)) {
    output$representator <- representotar
    output$distance <- distance
    output$cluster <- cluster
    output$rf <- list(rf)
    output$other <- list(other)
  } else {
    output$representator <- c(output$representator, representotar)
    output$distance <- c(output$distance, distance)
    output$cluster <- c(output$cluster, cluster)
    output$rf <- append(output$accuracy, list(rf))
    output$other <- append(output$other, list(other))
  }
  output
}

output_pre <- function(output) {
  as_tibble(output) %>% 
    tidyr::unnest_longer(accuracy, values_to = "accuracy", indices_to = "accuracy_method") %>% 
    tidyr::unnest_longer(accuracy, values_to = "accuracy", indices_to = "rf_method") %>% 
    tidyr::unnest_wider(accuracy)
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
           "average" = "Combination"
    ),
    "Two-level"
  )
  paste0(representor, distance, cluster)
}



mcb_ <- function(orig, type, plot_type, target,
                 file_path, mar = NULL, conf.level) {
  orig_ <- orig %>%
    select(method, rmsse, batch) %>%
    tidyr::nest(rmsse = -c("method")) %>%
    mutate_at("rmsse", purrr::map, function(g) {
      g <- g %>% arrange(batch)
      if (type == "hierarchy") {
        return(list(sapply(g$rmsse, mean)))
      } else {
        return(list(do.call(c, g$rmsse)))
      }
    }) %>%
    tidyr::unnest(rmsse)
  
  rmsse <- do.call(cbind, orig_$rmsse)
  colnames(rmsse) <- orig_$method
  pdf(file_path, width=8, height = 6)
  if (plot_type == "custom") {
    stopifnot(!is.null(mar))
    par(mar = mar)
    nemenyi(rmsse, plottype = "vmcb", target = target, conf.level = conf.level)
  }
  if (plot_type == "orig") {
    tsutils::nemenyi(rmsse, plottype = "vmcb", conf.level = conf.level)
  }
  dev.off()
}
mcb_series_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>%
    arrange(batch) %>%
    pull(rmsse) %>%
    do.call(c, .)
  rand_rmsse <-
    rand %>%
    arrange(batch, cluster) %>%
    select(cluster, rmsse, batch) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    mutate_at("rmsse", purrr::map, function(g) {
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
