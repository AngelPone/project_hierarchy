batch <- 0
path <- "mortality"
bfmethod <- "arima"

source("R/construct_hierarchy.R", chdir = T)

num.cores <- 6
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

if (!file.exists(store_path)) {
  stop("run configuration and base forecast first!")
}

BASEFORECAST_STORE <- readRDS(store_path)$bfstore
output <- readRDS(store_path)$output
data <- readRDS(store_path)$data

# natural hierarchy
S <- read.csv("data/S.csv", row.names = 1)
colnames(S) <- stringi::stri_replace_all_fixed(colnames(S), ".", "-")
natural_data <- data
natural_data$nl <- list(list(S = unname(as.matrix(S[,colnames(mortality)][which(!(rownames(S) %in% colnames(S))), ])),
                             basef = NULL,
                             rf = NULL))

natural_data <- forecast(natural_data, f.arima)
natural_data <- reconcile.all(natural_data)
accs <- evaluate.hts(natural_data, metrics, type = "nl")
output <- add_result(output, "", "", "natural", accs, other = list(
  S = as(natural_data$nl[[1]]$S, "sparseMatrix")
))



# random hierarchy
for (n_cluster in 3:10) {
  data <- build_level(data, representator.ts, 
              distance.euclidean, cluster.random, 
              n_clusters = rep(n_cluster, 20))
  data <- forecast(data, f.arima)
  data <- reconcile.all(data)
  accs <- evaluate.hts(data, metrics, type = "average")
  accs_single <- evaluate.hts(data, metrics, type = "nl")
  for (acc in seq_along(accs_single)) {
    output <- add_result(output, "", "", "random-single", accs_single[[acc]],
                         other = list(n_clusters = n_cluster, S = as(accs$nl[[acc]]$S, "sparseMatrix")))
  }
  
  output <- add_result(output, "", "", "random", accs, other = list(n_clusters = n_cluster, average=20,
                                                                    S = lapply(data$nl, function(x){ as(x$S, "sparseMatrix") })))
}

saveResult()







