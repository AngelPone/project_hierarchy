rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
batch <- as.integer(args[[1]])
path <- "tourism"
bfmethod <- "ets"
set.seed(810975)
source("R/construct_hierarchy.R", chdir = T)

num.cores <- 8
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

store_path <- sprintf("%s/store_%s.rds", path, batch)

BASEFORECAST_STORE <- readRDS(store_path)$bfstore
output <- readRDS(store_path)$output
data <- readRDS(store_path)$data
FEATURES <- readRDS(store_path)$features
DISTANCEMAT <- readRDS(store_path)$distance

S <- read.csv("tourism/S.csv", row.names = 1)
S <- unname(as.matrix(S))

# random hierarchy

for (n_repeat in c(100)) {
  for (n_cluster in 5:15) {
    print(sprintf("%s random hierarchy: %s repeats * %s cluster", as.character(Sys.time()), n_repeat, n_cluster))
    data <- build_level(data, "ts", 
                        "euclidean", cluster.random, 
                        n_clusters = rep(n_cluster, n_repeat))
    data <- forecast(data, bfmethod, frequency = 12, h=12)
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

saveResult()


# random natural
print(paste0(Sys.time(), "random natural 100 ......."))
new_S <- S[2:251,]
data$nl <- lapply(1:100, function(i) {
  list(S = new_S[,sample(NCOL(new_S), NCOL(new_S))],
       basef = NULL,
       rf = NULL)
})

data <- forecast(data, bfmethod, h=12, frequency = 12)
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







