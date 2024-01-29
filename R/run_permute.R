args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

set.seed(20231019)


source("R/construct_hierarchy.R", chdir = T)


dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- 12
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

for (batch in 0:(batch_length-1)) {
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
  
  data <- readRDS(store_path)
  tbl <- nl2tibble(data$nl)
  natural_S <- tbl$S[which(tbl$cluster=="natural")]
  stopifnot(length(natural_S) == 1)
  natural_S <- natural_S[[1]]
  
  # permute natural
  for (i in 1:100){
    new_S <- natural_S[,sample(NCOL(natural_S))]
    data <- add_nl(data, new_S, "", "", paste0("permute-natural-", i))
  }
  
  # permute clustering
  
  for (distance in DISTANCES) {
    for (representor in REPRESENTORS) {
      for (cluster in c("Kmedoids-dr", "hcluster-dr")) {
        S <- tbl$S[which((tbl$representor==representor) & (tbl$cluster == cluster) & (tbl$representor==representor))]
        stopifnot(length(S) == 1)
        S <- S[[1]]
        for (i in 1:100) {
          new_S <- S[,sample(NCOL(S))]
          data <- add_nl(data, new_S, representor, distance, paste0("permute-",cluster,"-", i))
        }
      }
    }
  }
  
  for (distance in DISTANCES2) {
    for (representor in REPRESENTORS2) {
      for (cluster in c("Kmedoids", "hcluster-dr")) {
        S <- tbl$S[which((tbl$representor==representor) & (tbl$cluster == cluster) & (tbl$representor==representor))]
        stopifnot(length(S) == 1)
        S <- S[[1]]
        for (i in 1:100) {
          new_S <- S[,sample(NCOL(S))]
          data <- add_nl(data, new_S, representor, distance, paste0("permute-",cluster,"-", i))
        }
      }
    }
  }
  saveRDS(data, store_path)
}

