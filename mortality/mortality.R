path <- "mortality"
bfmethod <- "ets"
num.cores <- 8
set.seed(42)

cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

source("R/construct_hierarchy.R", chdir = T)


mortality <- read.csv("data/mortality.csv")
colnames(mortality) <- stringi::stri_replace_all_fixed(colnames(mortality), ".", "-")


for (batch in 0:11) {
  
  store_path <- sprintf("mortality/store_%s.rds", batch)
  
  data <- hts(rbind(rep(1, 98), diag(98)), 
              bts = unname(as.matrix(mortality))[1:(240-batch),],
              tts = unname(as.matrix(mortality))[(241 - batch):(252 - batch),])
  # BASEFORECAST_STORE <- list()
  # BASEFORECAST_STORE[[bfmethod]] <- list()
  BASEFORECAST_STORE <- readRDS(sprintf("mortality/ets_backup/store_%s.rds", batch))$bfstore
  output <- create_new_output()
  
  ## bottom and total forecast
  data <- forecast(data, bfmethod, h=12, frequency=12)
  accs <- evaluate.hts(data, metrics, type = "base")
  output <- add_result(output, "", "", "", accs)
  data <- reconcile.all(data)
  accs <- evaluate.hts(data, metrics, type = "rf")
  output <- add_result(output, "", "", "", accs)
  
  
  print("computing features .....")
  # FEATURES <- NULL
  FEATURES <- readRDS(sprintf("mortality/ets_backup/store_%s.rds", batch))$features
  # features.compute(data, frequency = 12)
  
  print("computing distance matrix ....")
  # DISTANCEMAT <- list()
  # for (representor in c("ts", "error", "ts.features", "error.features", "forecast")) {
  #   DISTANCEMAT[[representor]] <- list()
  #   for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
  #     cluster_input <- get(paste0("representator.", representor))(data)
  #     distance_mat <- matrix(0, 98, 98)
  #     
  #     for(i in 1:98) {
  #       for (j in 1:i) {
  #         distance_mat[i, j] <- get(paste0("distance.", distance))(cluster_input[, i], cluster_input[, j])
  #         distance_mat[j, i] <- distance_mat[i, j]
  #       }
  #     }
  #     DISTANCEMAT[[representor]][[distance]] <- distance_mat
  #   }
  # }
  DISTANCEMAT <- readRDS(sprintf("mortality/ets_backup/store_%s.rds", batch))$distance
  
  saveResult()
}

