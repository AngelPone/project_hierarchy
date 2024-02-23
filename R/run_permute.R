args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- "ets"

set.seed(20231019)


source("R/construct_hierarchy.R", chdir = T)

dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- as.integer(args[[2]])
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

REPRESENTORS <- c("ts-dr", "error-dr", "ts.features-dr", "error.features-dr")
DISTANCES <- rep("euclidean", 4)

REPRESENTORS2 <- c("ts", "error")
DISTANCES2 <- rep("dtw", 2)

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

get_best_clustering <- function(path, batch_length, h) {
  all_rmsse <- NULL
  for (batch in 0:(batch_length - 1)) {
    store_path <- sprintf("%s/ets/batch_%s.rds", path, batch)
    data <- readRDS(store_path)
    cluster <- sapply(data$nl, function(x) {
      x$cluster
    })
    idx <- which((!startsWith(cluster, "permute") & (cluster != "natural") & (cluster != "")))
    stopifnot(length(idx) == 12)
    tts <- cbind(rowSums(data$tts), data$tts)
    bts <- cbind(rowSums(data$bts), data$bts)
    tbl <- nl2tibble(data$nl[idx])
    tbl$rmsse <- sapply(tbl$rf, function(x) {
      mean(sapply(1:NCOL(x$mint), function(f) {
        metric.rmsse(tts[1:h, f], x$mint[1:h, f], bts[, f])
      }))
    })
    tbl <- select(tbl, representor, cluster, distance, rmsse) %>%
      mutate(batch = batch)
    all_rmsse <- rbind(all_rmsse, tbl)
  }
  all_rmsse %>%
    group_by(representor, distance, cluster) %>%
    summarise(rmsse = mean(rmsse)) %>%
    arrange(rmsse)
}
print("Calculating best performing method ....")
performing_sort <- get_best_clustering(path,
  batch_length,
  h = forecast_horizon
)
print("best performing method:")
print(performing_sort[1, ])

new_permute <- lapply(1:100, function(x) {
  sample(m)
})


for (batch in 0:batch_length) {
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)

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
    representor <- performing_sort$representor[[1]]
    distance <- performing_sort$distance[[1]]
    cluster <- performing_sort$cluster[[1]]
    data <- generate_randomization(data, representor, distance, cluster, tbl, new_permute)
    stopifnot((length(data$nl) - 14) %% 100 == 0)
  }
  saveRDS(data, store_path)
}
