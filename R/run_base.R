args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

num.cores <- 8

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
  
  data <- hts(rbind(rep(1, m), diag(m)),
               bts = dt$data[1:(96+batch), (n-m+1):n],
               tts = dt$data[(96+batch+1):(96+batch+forecast_horizon), (n-m+1):n])
  
  data <- hts.basef(data, bfmethod, h=forecast_horizon, frequency=12)
  
  print(paste0(Sys.time(), "computing features ..."))
  data <- features.compute(data, frequency = frequency)
  
  print(paste0(Sys.time(), "computing distance matrix ..."))
  DISTANCEMAT <- list()  
  for (representor in c("ts", "error", "ts.features", "error.features", "forecast")) {
    DISTANCEMAT[[representor]] <- list()
    for (distance in c("euclidean", "dtw", "negcor", "cor", "uncorrelation")) {
      cluster_input <- get(paste0("representator.", representor))(data)
      distance_method <- get(paste0("distance.", distance))
      distance_mat <- matrix(0, m, m)

      lst <- foreach(row=1:m, .packages = c("dtw")) %dopar% {
        output <- c()
        for (col in 1:row) {
          dis <- distance_method(cluster_input[, row], cluster_input[, col])
          output <- c(output, dis)
        }
        output
      }
      for (row in 1:m) {
        distance_mat[row, 1:row] <- lst[[row]]
        distance_mat[1:row, row] <- lst[[row]]
      }
      
      DISTANCEMAT[[representor]][[distance]] <- distance_mat
    }
  }
  data$distance <- DISTANCEMAT
  
  data <- add_nl(data, NULL, "", "", "")
  data <- add_nl(data, dt$S[2:(n-m),], representor = "", distance = "", cluster = "natural")
  saveRDS(data, store_path)
}




