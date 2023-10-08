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


REPRESENTORS <- c("ts", "error", "ts.features", "error.features", "forecast")
DISTANCES <- c("euclidean", "dtw", "negcor", "cor", "uncorrelation")


# load dataset
for (batch in 0:batch_length) {
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
  data <- readRDS(store_path)
  
  print(paste0(Sys.time(), " KMedoids ..."))
  for (representor in REPRESENTORS) {
    for (distance in DISTANCES) {
      nl <- foreach(n_clusters = 5:15, .packages = "cluster") %dopar% {
        build_level(hts = data, representor=representor, 
                    distance = distance, 
                    cluster = cluster.kmedoids,
                    n_clusters = n_clusters)[[1]]
      }
      
      for (l in seq_along(nl)) {
        data <- add_nl(data, nl[[l]], representor, distance, paste0("Kmedoids-", l + 4))
      }
      
      unnested <- do.call(rbind,  nl)
      data <- add_nl(data, unnested, representor, distance, paste0("unnestedKmedoids"))
    }
  }
  
  print(paste0(Sys.time(), " hierarchical clustering ..."))
  
  for (representor in REPRESENTORS) {
    for (distance in DISTANCES) {
      for (method in c("ward", "average", "single", "complete", "weighted")) {
        nl <- build_level(hts = data, representor=representor, 
                          distance = distance, 
                          cluster = cluster.hcluster,
                          method = method)[[1]]
        data <- add_nl(data, nl, representor, distance, paste0("hcluster-", method))
      }
    }
  }
  
  print(paste0(Sys.time(), " nested Kmedoids 2 ..."))
  REPRESENTORS_comb <- rep(REPRESENTORS, length(DISTANCES))
  DISTANCES_comb <- rep(DISTANCES, each=length(REPRESENTORS))
  
  nlevel <- which(m / 2**c(3:6) < 10 & m / 2**c(3:6) > 5) + 2
  nl <- foreach(representor = REPRESENTORS_comb, distance = DISTANCES_comb,
                .packages = "cluster") %dopar% {
                  build_level(hts = data, representor=representor, 
                              distance = distance, 
                              cluster = cluster.nestedkmedoids,
                              n_clusters = rep(2, nlevel))
                }
  for (l in seq_along(nl)) {
    data <- add_nl(data, nl[[l]][[1]], REPRESENTORS_comb[l], DISTANCES_comb[l], 
                   paste0("nestedkmedoids-", do.call(paste0, as.list(rep(2, nlevel)))))
  }
  
  print(paste0(Sys.time(), " nested Kmedoids 3 ..."))
  
  nlevel <- which(m / 3**c(2:5) < 4 & m / 3**c(2:5) > 3) + 1
  nl <- foreach(representor = REPRESENTORS_comb, distance = DISTANCES_comb,
               .packages = "cluster") %dopar% {
            build_level(hts = data, representor=representor, 
                        distance = distance, 
                        cluster = cluster.nestedkmedoids,
                        n_clusters = rep(3, nlevel))
               }
  for (l in seq_along(nl)) {
    data <- add_nl(data, nl[[l]][[1]], REPRESENTORS_comb[l], DISTANCES_comb[l], 
                   paste0("nestedkmedoids-", do.call(paste0, as.list(rep(3, nlevel)))))
  }
  
  saveRDS(data, store_path)
}










