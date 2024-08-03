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

for (batch in 0:(batch_length-1)) {
  print(sprintf("%s batch %s ....", Sys.time(), batch))
  store_path <- sprintf("%s/batch_%s.rds", path, batch)
  data <- readRDS(store_path)

  data <- hts.nlf(data, frequency = frequency, h = forecast_horizon)

  saveRDS(data, store_path)
}

