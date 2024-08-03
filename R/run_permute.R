args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]

set.seed(20231019)


source("R/utils.R", chdir = T)

dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon + 1

REPRESENTORS <- c("ts-dr", "error-dr", "ts.features-dr", "error.features-dr")
DISTANCES <- rep("euclidean", 4)

REPRESENTORS2 <- c("ts", "error")
DISTANCES2 <- rep("dtw", 2)


generate_randomization <- function(data, representor, distance, cluster, tbl, permute) {
  S <- tbl$S[which((tbl$representor == representor) & (tbl$cluster == cluster) & (tbl$distance == distance))]
  permute_S <- tbl$S[which((tbl$representor == representor) & (startsWith(tbl$cluster, paste0("permute-", cluster, "-"))) & (tbl$distance == distance))]
  if (length(permute_S) == 100) {
    print("existing permutation, skip ...")
    return(data)
  }
  stopifnot(length(permute_S) == 0)
  stopifnot(length(S) == 1)
  S <- S[[1]]
  for (i in 1:100) {
    new_S <- S[, permute[[i]]]
    data <- add_nl(data, new_S, representor, distance, paste0("permute-", cluster, "-", i))
  }
  data
}


new_permute <- lapply(1:100, function(x) {
  sample(m)
})


for (batch in 0:(batch_length-1)) {
  store_path <- sprintf("%s/batch_%s.rds", path, batch)

  data <- readRDS(store_path)
  tbl <- nl2tibble(data$nl)

  natural_S <- tbl$S[[which(tbl$cluster == "natural")]]
  if (sum(startsWith(tbl$cluster, "permute-natural")) == 0) {
    for (i in 1:100) {
      data <- add_nl(data, natural_S[, new_permute[[i]]], "", "", paste0("permute-natural-", i))
    }
  }

  # permute clustering
  if (path == "mortality") {
    for (i in seq_along(DISTANCES)) {
      representor <- REPRESENTORS[i]
      distance <- DISTANCES[i]
      for (cluster in c("Kmedoids-dr", "hcluster-dr")) {
        data <- generate_randomization(data, representor, distance, cluster, tbl, new_permute)
      }
    }

    for (i in seq_along(DISTANCES2)) {
      representor <- REPRESENTORS2[i]
      distance <- DISTANCES2[i]
      for (cluster in c("Kmedoids", "hcluster")) {
        data <- generate_randomization(data, representor, distance, cluster, tbl, new_permute)
      }
    }
    stopifnot(length(data$nl) == 1314)
  }
  if (path == "tourism") {
    best_ <- readRDS("tourism/eval_cluster.rds")$best_
    data <- generate_randomization(data, best_$representor, best_$distance, best_$cluster, tbl, new_permute)
    stopifnot(length(data$nl) == 214)
  }
  
  
  saveRDS(data, store_path)
}
