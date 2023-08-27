rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
batch <- as.integer(args[[1]])
path <- "mortality"
bfmethod <- "arima"
set.seed(42)

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

# natural hierarchy
S <- read.csv("data/S.csv", row.names = 1)
colnames(S) <- stringi::stri_replace_all_fixed(colnames(S), ".", "-")
natural_data <- data

natural_data$nl <- list(list(S = unname(as.matrix(S[which(!(rownames(S) %in% colnames(S))), ])),
                             basef = NULL,
                             rf = NULL))

natural_data <- forecast(natural_data, bfmethod)
natural_data <- reconcile.all(natural_data)
accs <- evaluate.hts(natural_data, metrics, type = "nl")
output <- add_result(output, "", "", "natural", accs, other = list(
  S = as(natural_data$nl[[1]]$S, "sparseMatrix")
))



# random hierarchy

for (n_repeat in c(20, 50, 100)) {
  for (n_cluster in 3:10) {
    data <- build_level(data, representator.ts, 
                        distance.euclidean, cluster.random, 
                        n_clusters = rep(n_cluster, n_repeat))
    data <- forecast(data, bfmethod)
    data <- reconcile.all(data)
    accs <- evaluate.hts(data, metrics, type = "average")
    accs_single <- evaluate.hts(data, metrics, type = "nl")
    
    accs_single <- lapply(1:n_repeat, function(x) {
      tmp <- lapply(metrics, function(metric){
        accs_single[[metric]][[x]]
      })
      names(tmp) <- metrics
      tmp
    })
    
    for (acc in seq_along(accs_single)) {
      output <- add_result(output, "", "", paste0("random-single-", n_cluster), accs_single[[acc]],
                           other = list(n_clusters = n_cluster, S = as(data$nl[[acc]]$S, "sparseMatrix")))
    }
    
    output <- add_result(output, "", "", paste0("random-average-", n_cluster, "-", n_repeat), accs, other = list(n_clusters = n_cluster, average=n_repeat,
                                                                                                                 S = lapply(data$nl, function(x){ as(x$S, "sparseMatrix") })))
  }
}


# random natural
new_S <- unname(as.matrix(S[which(!(rownames(S) %in% colnames(S))), ]))
data$nl <- lapply(1:20, function(i) {
  list(S = new_S[,sample(NCOL(new_S), NCOL(new_S))],
            basef = NULL,
            rf = NULL)
})

data <- forecast(data, "arima")
data <- reconcile.all(data)
accs <- evaluate.hts(data, metrics, type="average")
accs_single <- evaluate.hts(data, metrics, type="nl")

accs_single <- lapply(1:20, function(x) {
  tmp <- lapply(metrics, function(metric){
    accs_single[[metric]][[x]]
  })
  names(tmp) <- metrics
  tmp
})

for (acc in seq_along(accs_single)) {
  output <- add_result(output, "", "", paste0("random-single-natural"), accs_single[[acc]],
                       other = list(S = as(data$nl[[acc]]$S, "sparseMatrix")))
}

output <- add_result(output, "", "", paste0("random-average-natural-", 20), accs, other = list(average=20, S = lapply(data$nl, function(x){ as(x$S, "sparseMatrix") })))


data$nl <- lapply(1:100, function(i) {
  list(S = new_S[,sample(NCOL(new_S), NCOL(new_S))],
       basef = NULL,
       rf = NULL)
})

data <- forecast(data, "arima")
data <- reconcile.all(data)
accs <- evaluate.hts(data, metrics, type="average")
accs_single <- evaluate.hts(data, metrics, type="nl")

accs_single <- lapply(1:100, function(x) {
  tmp <- lapply(metrics, function(metric){
    accs_single[[metric]][[x]]
  })
  names(tmp) <- metrics
  tmp
})

for (acc in seq_along(accs_single)) {
  output <- add_result(output, "", "", paste0("random-single-natural"), accs_single[[acc]],
                       other = list(S = as(data$nl[[acc]]$S, "sparseMatrix")))
}

output <- add_result(output, "", "", paste0("random-average-natural-", 100), accs, other = list(average=100, S = lapply(data$nl, function(x){ as(x$S, "sparseMatrix") })))

saveResult()







