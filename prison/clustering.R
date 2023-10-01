
path <- "prison"
bfmethod <- "arima"
batch <- 1
source("R/construct_hierarchy.R", chdir = T)
set.seed(42)
# rm(.Random.seed, envir=globalenv())
num.cores <- 8
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)


for (batch in 0:7) {
  
  store_path <- sprintf("prison/store_%s.rds", batch)
  
  BASEFORECAST_STORE <- readRDS(store_path)$bfstore
  output <- readRDS(store_path)$output
  data <- readRDS(store_path)$data
  FEATURES <- readRDS(store_path)$features
  
  if (is.null(FEATURES)) {
    features.compute(data, frequency = 4)
  }
  
  # K-medoids
  for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
    for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
      print(sprintf("%s-%s", representator, distance))
      for (n_clusters in 2:5) {
        dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                          distance = get(paste0("distance.", distance)), 
                          cluster = cluster.kmedoids,
                          n_clusters = n_clusters)
        dt <- forecast(dt, bfmethod, h=NROW(data$tts), frequency = 4)
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
        dt <- build_level(hts = data, representator = get(paste0("representator.", representator)),
                          distance = get(paste0("distance.", distance)),
                          cluster = cluster.hcluster,
                          method = method)
        dt <- forecast(dt, bfmethod, h=NROW(data$tts), frequency = 4)
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
      print(sprintf("nested 22 %s-%s", representator, distance))
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                        distance = get(paste0("distance.", distance)), 
                        cluster = cluster.nestedkmedoids,
                        n_clusters = c(2, 2, 2))
      dt <- forecast(dt, bfmethod, h=NROW(data$tts), frequency = 4)
      dt <- reconcile.all(dt)
      accs <- evaluate.hts(dt, metrics, type = "nl")
      output <- add_result(output, 
                           representotar = representator, 
                           distance = distance, 
                           cluster = "kmedoids-d-nested-222", 
                           accs = accs,
                           other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
    }
  }
  
  saveResult()
  
  for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
    for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
      print(sprintf("nested 333 %s-%s", representator, distance))
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)), 
                        distance = get(paste0("distance.", distance)), 
                        cluster = cluster.nestedkmedoids,
                        n_clusters = c(3, 3))
      dt <- forecast(dt, bfmethod, h=NROW(data$tts), frequency = 4)
      dt <- reconcile.all(dt)
      accs <- evaluate.hts(dt, metrics, type = "nl")
      output <- add_result(output, 
                           representotar = representator, 
                           distance = distance, 
                           cluster = "kmedoids-d-nested-33", 
                           accs = accs,
                           other = list(S = as(dt$nl[[1]]$S, "sparseMatrix")))
    }
  }
  
  saveResult()
  
  
  # unnested
  for (representator in c("ts", "error", "ts.features", "error.features", "forecast")) {
    for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
      print(sprintf("unnested %s-%s", representator, distance))
      dt <- build_level(hts = data, representator = get(paste0("representator.", representator)),
                        distance = get(paste0("distance.", distance)),
                        cluster = cluster.kmedoids,
                        n_clusters = 2:5)
      dt <- forecast(dt, bfmethod, h=NROW(data$tts), frequency = 4)
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
}

