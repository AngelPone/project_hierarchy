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



