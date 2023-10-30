args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

library(cluster)
library(foreach)
source("R/representator.R")

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

cal_gap <- function(dt, representor) {
  a <- clusGap(t(representor(dt, dr = dr)), pam, K.max = 50, B=50)
  maxSE(a$Tab[,3], a$Tab[,4], method = "Tibs2001SEmax")
}



dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

REPRESENTORS <- paste0(c("ts", "error", "ts.features", "error.features"), "-dr")

for (batch in 0:batch_length) {
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
  data <- readRDS(store_path)
  nl <- nl2tibble(data$nl)
  print(sprintf("%s %s", Sys.time(), batch))
  
  drts <- list()
  for (representor in REPRESENTORS) {
    rsplit <- strsplit(representor, "-")[[1]]
    dr <- rsplit[2] == "dr"
    drts[[representor]] <- t(get(paste0("representator.", rsplit[1]))(data, dr=dr))
  }

  drtgaps <- foreach(drt = iterators::iter(drts), .packages = "cluster") %dopar% {
    a <- clusGap(drt, pam, K.max = 50, B=50)
    maxSE(a$Tab[,3], a$Tab[,4], method = "Tibs2001SEmax")
  }
  names(drtgaps) <- names(drts)
  
  for (representor in REPRESENTORS) {
    idx <- which(nl$representor == representor & nl$cluster == "Kmedoids-dr")
    data$nl[[idx]]$other$gap <- drtgaps[[representor]]
  }
  
  saveRDS(data, store_path)
}