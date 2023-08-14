rm(list = ls())
source("R/construct_hierarchy.R", chdir = T)
library(Matrix)

num.cores <- 8
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

data <- readRDS("mortality/data.rds")


# K-medoids
for (representator in c("ts", "error")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor")) {
    for (n_clusters in 3:10) {
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                        distance = get(paste0("distance.", distance)), 
                        cluster = cluster.kmedoids,
                        n_clusters = n_clusters)
      dt <- forecast(dt, f.arima)
      dt <- reconcile.all(dt)
      accs <- evaluate.hts(dt, metrics, type = "nl")
      output <- add_result(output, 
                           representotar = representator, 
                           distance = distance, 
                           cluster = "kmedoids", 
                           accs = accs,
                           other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
    }
  }
}


for (representator in c("ts", "error")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor")) {
    for (method in c("ward", "average", "single", "complete", "weighted")) {
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                        distance = get(paste0("distance.", distance)), 
                        cluster = cluster.hcluster,
                        method = method)
      dt <- forecast(dt, f.arima)
      dt <- reconcile.all(dt)
      accs <- evaluate.hts(dt, metrics, type = "nl")
      output <- add_result(output, 
                           representotar = representator, 
                           distance = distance, 
                           cluster = paste0("hcluster-", method), 
                           accs = accs,
                           other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
    }
  }
}

saveRDS(output, "mortality/ts_error.rds")








