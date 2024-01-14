args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

num.cores <- 8
source("R/construct_hierarchy.R", chdir = T)
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)


dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

for (batch in 0:(batch_length-1)) {
  print(sprintf("%s batch %s ....", Sys.time(), batch))
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
  data <- readRDS(store_path)
  
  data <- hts.nlf(data, bfmethod, frequency=frequency, h=forecast_horizon)
  
  saveRDS(data, store_path)
}