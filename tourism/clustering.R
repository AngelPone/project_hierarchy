args <- commandArgs(trailingOnly = TRUE)
batch <- as.integer(args[[1]])
path <- "tourism"
bfmethod <- "ets"
source("R/construct_hierarchy.R", chdir = T)
set.seed(42)
# rm(.Random.seed, envir=globalenv())
num.cores <- 8
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

store_path <- sprintf("tourism/store_%s.rds", batch)

BASEFORECAST_STORE <- readRDS(store_path)$bfstore
output <- readRDS(store_path)$output
data <- readRDS(store_path)$data
FEATURES <- readRDS(store_path)$features
DISTANCEMAT <- readRDS(store_path)$distance

# K-medoids
for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
    print(sprintf("Kmedoids: %s-%s", representator, distance))
    for (n_clusters in 5:15) {
      dt <- build_level(hts = data, representator = representator,
                        distance = distance,
                        cluster = cluster.kmedoids,
                        n_clusters = n_clusters)
      dt <- forecast(dt, bfmethod, frequency = 12, h = 12)
      dt <- reconcile.all(dt)
      accs <- evaluate.hts(dt, metrics, type = "nl")
      output <- add_result(output,
                           representotar = representator,
                           distance = distance,
                           cluster = paste0("kmedoids-", n_clusters),
                           accs = accs,
                           other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
    }
  }
}

saveResult()

# hierarchical clustering
for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
    for (method in c("ward", "average", "single", "complete", "weighted")) {
      print(sprintf("%s, %s", distance, method))
      dt <- build_level(hts = data, representator = representator,
                        distance = distance,
                        cluster = cluster.hcluster,
                        method = method)
      dt <- forecast(dt, bfmethod, frequency = 12, h = 12)
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
saveResult()

# K-medoids, nested
for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
    print(sprintf("nested 2222 %s-%s", representator, distance))
    dt <- build_level(hts = data, representator = representator, 
                      distance = distance, 
                      cluster = cluster.nestedkmedoids,
                      n_clusters = c(2, 2, 2, 2,2))
    dt <- forecast(dt, bfmethod, frequency = 12, h = 12)
    dt <- reconcile.all(dt)
    accs <- evaluate.hts(dt, metrics, type = "nl")
    output <- add_result(output, 
                         representotar = representator, 
                         distance = distance, 
                         cluster = "kmedoids-d-nested-22222", 
                         accs = accs,
                         other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
  }
}

saveResult()

for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
    print(sprintf("nested 333 %s-%s", representator, distance))
    dt <- build_level(hts = data, representator = representator, 
                      distance = distance, 
                      cluster = cluster.nestedkmedoids,
                      n_clusters = c(3, 3, 3, 3))
    dt <- forecast(dt, bfmethod, frequency = 12, h = 12)
    dt <- reconcile.all(dt)
    accs <- evaluate.hts(dt, metrics, type = "nl")
    output <- add_result(output, 
                         representotar = representator, 
                         distance = distance, 
                         cluster = "kmedoids-d-nested-3333", 
                         accs = accs,
                         other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
  }
}

saveResult()


# unnested
for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
  for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
    print(sprintf("unnested %s-%s", representator, distance))
      dt <- build_level(hts = data, representator = representator,
                        distance = distance,
                        cluster = cluster.kmedoids,
                        n_clusters = 3:20)
      dt <- forecast(dt, bfmethod, frequency = 12, h = 12)
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


