args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]

source("R/utils.R", chdir = T)


dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon + 1


REPRESENTORS <- c("ts-dr", "error-dr", "ts.features-dr", "error.features-dr", "ts", "error")
DISTANCES <- c(rep("euclidean", 4), "dtw", "dtw")


for (batch in 0:(batch_length - 1)) {
  # path to save batch
  store_path <- sprintf("%s/batch_%s.rds", path, batch)
  
  # construct hts
  data <- hts(rbind(rep(1, m), diag(m)),
              bts = dt$data[1:(96 + batch), (n - m + 1):n],
              tts = dt$data[(96 + batch + 1):(96 + batch + forecast_horizon), (n - m + 1):n]
  )

  # produce base forecast
  data <- hts.basef(data, h = forecast_horizon, frequency = 12)

  # compute features
  print(paste0(Sys.time(), " computing features ..."))
  data <- features.compute(data, frequency = frequency)

  # compute distance matrix for clustering
  print(paste0(Sys.time(), "computing distance matrix ..."))
  DISTANCEMAT <- list()
  for (i in seq_along(REPRESENTORS)) {
    representor <- REPRESENTORS[i]
    distance <- DISTANCES[i]
    cluster_input <- get(paste0("representor.", strsplit(representor, "-")[[1]][1]))(data, dr = endsWith(representor, "-dr"))
    distance_method <- get(paste0("distance.", distance))
    distance_mat <- matrix(0, m, m)
    lst <- future_map(1:m, function(row) {
      output <- c()
      for (col in 1:row) {
        dis <- distance_method(cluster_input[, row], cluster_input[, col])
        output <- c(output, dis)
      }
      output
    })
    for (row in 1:m) {
      distance_mat[row, 1:row] <- lst[[row]]
      distance_mat[1:row, row] <- lst[[row]]
    }
    DISTANCEMAT[[representor]][[distance]] <- distance_mat
    data$distance <- DISTANCEMAT
  }
  
  data <- add_nl(data, NULL, "", "", "")
  data <- add_nl(data, dt$S[2:(n - m), ], representor = "", distance = "", cluster = "natural")
  saveRDS(data, store_path)
}
