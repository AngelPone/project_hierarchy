
#' Euclidean distance
distance.euclidean <- function(x, y) { sqrt(sum((x - y)^2)) }


library(dtw)
#' Dynamic Time Warping
distance.dtw <- function(x, y) { dtw(x, y, distance.only = TRUE)$distance }


#' Correlation
distance.negcor <- function(x, y) { cov(x, y) + 1 }
distance.cor <- function(x, y) { abs(cov(x, y) - 1) }
