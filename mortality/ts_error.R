rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
batch <- as.integer(args[[1]])
path <- "mortality"
bfmethod <- "arima"
source("R/construct_hierarchy.R", chdir = T)

num.cores <- 8
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

if (!file.exists(store_path)) {
  stop("run configuration and base forecast first!")
}

BASEFORECAST_STORE <- readRDS(store_path)$bfstore
output <- readRDS(store_path)$output
data <- readRDS(store_path)$data
FEATURES <- readRDS(store_path)$features



# K-medoids

dt2 <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                   distance = get(paste0("distance.", distance)), 
                   cluster = cluster.kmedoids,
                   n_clusters = 6)

for (representator in c("ts", "error")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor")) {
    for (n_clusters in 3:10) {
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                        distance = get(paste0("distance.", distance)), 
                        cluster = cluster.kmedoids,
                        n_clusters = n_clusters)
      dt <- forecast(dt, bfmethod)
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
      dt <- forecast(dt, bfmethod)
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


representator <- "ts"
distance <- "euclidean"

# K-medoids, nested
for (representator in c("ts", "error")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor")) {
    dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                      distance = get(paste0("distance.", distance)), 
                      cluster = cluster.nestedkmedoids,
                      n_clusters = c(2, 2, 2, 2))
    dt <- forecast(dt, bfmethod)
    dt <- reconcile.all(dt)
    accs <- evaluate.hts(dt, metrics, type = "nl")
    output <- add_result(output, 
                         representotar = representator, 
                         distance = distance, 
                         cluster = "kmedoids-d-nested", 
                         accs = accs,
                         other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
  }
}

for (representator in c("ts", "error")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor")) {
    dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                      distance = get(paste0("distance.", distance)), 
                      cluster = cluster.nestedkmedoids,
                      n_clusters = c(3, 3, 3))
    dt <- forecast(dt, bfmethod)
    dt <- reconcile.all(dt)
    accs <- evaluate.hts(dt, metrics, type = "nl")
    output <- add_result(output, 
                         representotar = representator, 
                         distance = distance, 
                         cluster = "kmedoids-d-nested", 
                         accs = accs,
                         other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
  }
}



for (representator in c("ts", "error")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor")) {
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                        distance = get(paste0("distance.", distance)), 
                        cluster = cluster.kmedoids,
                        n_clusters = 3:20)
      dt <- forecast(dt, bfmethod)
      dt <- reconcile.all(dt)
      accs <- evaluate.hts(dt, metrics, type = "nl")
      output <- add_result(output, 
                           representotar = representator, 
                           distance = distance, 
                           cluster = "kmedoids-d-unnested", 
                           accs = accs,
                           other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
  }
}


saveResult()


