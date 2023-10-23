args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

num.cores <- 8
set.seed(20231019)

cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

source("R/construct_hierarchy.R", chdir = T)


dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

# load dataset
for (batch in 0:batch_length) {
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
  data <- readRDS(store_path)
  
  for (i in 1:50) {
    for (ncluster in 5:15) {
      S <- cluster.random(data$distance$ts$euclidean, n_clusters = ncluster)[[1]]
      data <- add_nl(data, S, "", "", paste0("random-", ncluster))
    }
  }
  
  newS <- dt$S[2:(n-m), ]
  
  for (i in 1:50) {
    S <- newS[, sample(m, m)]
    data <- add_nl(data, S, "", "", "random-natural")
  }
  
  saveRDS(data, store_path)
}










