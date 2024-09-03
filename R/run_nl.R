args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]

source("R/utils.R")

REPRESENTORS <- c("ts-dr", "error-dr", "ts.features-dr", "error.features-dr")
DISTANCES <- rep("euclidean", 4)


REPRESENTORS2 <- c("ts", "error")
DISTANCES2 <- rep("dtw", 2)


for (batch in 0:(batch_length-1)) {
  store_path <- sprintf("%s/batch_%s.rds", path, batch)
  data <- readRDS(store_path)

  for (i in seq_along(REPRESENTORS)) {
    nl <- build_level(
      hts = data, representor = REPRESENTORS[i],
      distance = DISTANCES[i],
      cluster = cluster.kmedoids,
      n_clusters = 1:(m - 1)
    )
    data <- add_nl(data, nl$S, REPRESENTORS[i], DISTANCES[i], "Kmedoids-dr",
      other = nl$info
    )
  }

  for (i in seq_along(REPRESENTORS2)) {
    nl <- build_level(
      hts = data, representor = REPRESENTORS2[i],
      distance = DISTANCES2[i],
      cluster = cluster.kmedoids,
      n_clusters = 1:(m - 1)
    )
    data <- add_nl(data, nl$S, REPRESENTORS2[i], DISTANCES2[i], "Kmedoids",
      other = nl$info
    )
  }

  # hierarchical clustering dimension reduction
  for (i in seq_along(REPRESENTORS)) {
    nl <- build_level(
      hts = data, representor = REPRESENTORS[i],
      distance = DISTANCES[i],
      cluster = cluster.hcluster,
      method = "ward"
    )[[1]]
    data <- add_nl(data, nl, REPRESENTORS[i], DISTANCES[i], "hcluster-dr")
  }

  for (i in seq_along(REPRESENTORS2)) {
    print(sprintf("%s hcluster %s * %s ", Sys.time(), REPRESENTORS2[i], DISTANCES2[i]))
    nl <- build_level(
      hts = data, representor = REPRESENTORS2[i],
      distance = DISTANCES2[i],
      cluster = cluster.hcluster,
      method = "ward"
    )[[1]]
    data <- add_nl(data, nl, REPRESENTORS2[i], DISTANCES2[i], "hcluster")
  }

  saveRDS(data, store_path)
}
