# configuration

path <- "tourism"
bfmethod <- "ets"
num.cores <- 8
set.seed(42)

cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

source("R/construct_hierarchy.R", chdir = T)


# dt <- readRDS("tourism/data.rds")

## load dataset
for (batch in 0:11) {
  store_path <- sprintf("tourism/store_%s.rds", batch)
  
  # data <- hts(rbind(rep(1, 304), diag(304)),
  #             bts = dt$data[1:(216-batch), 252:555],
  #             tts = dt$data[(217-batch):(228-batch), 252:555])
  # 
  # ## load base forecast store
  # BASEFORECAST_STORE <- list()
  # BASEFORECAST_STORE[[bfmethod]] <- list()
  orig_path <- sprintf("tourism/ets_backup/store_%s.rds", batch)
  
  BASEFORECAST_STORE <- readRDS(orig_path)$bfstore
  FEATURES <- readRDS(orig_path)$features
  DISTANCEMAT <- readRDS(orig_path)$distance
  
  output <- create_new_output()
  data <- readRDS(orig_path)$data
  data$basef <- NULL
  data$rf <- NULL
  data$resid <- NULL
  
  ## bottom and total forecast
  data <- forecast(data, bfmethod, h=12, frequency = 12)
  accs <- evaluate.hts(data, metrics, type = "base")
  output <- add_result(output, "", "", "", accs)
  data <- reconcile.all(data)
  accs <- evaluate.hts(data, metrics, type = "rf")
  output <- add_result(output, "", "", "", accs)
  # FEATURES <- NULL
  # features.compute(data, frequency = 12)
  
  
  # output <- readRDS(store_path)$output
  
 
  
  # distance matrix
  # print("computing distance matrix ....")
  # DISTANCEMAT <- list()
  # for (representor in c("ts", "error", "ts.features", "error.features", "forecast")) {
  #   DISTANCEMAT[[representor]] <- list()
  #   for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
  #     cluster_input <- get(paste0("representator.", representor))(data)
  #     distance_mat <- matrix(0, 304, 304)
  #     
  #     for(i in 1:304) {
  #       for (j in 1:i) {
  #         distance_mat[i, j] <- get(paste0("distance.", distance))(cluster_input[, i], cluster_input[, j])
  #         distance_mat[j, i] <- distance_mat[i, j]
  #       }
  #     }
  #     DISTANCEMAT[[representor]][[distance]] <- distance_mat
  #   }
  # }
  
  saveResult()
}




