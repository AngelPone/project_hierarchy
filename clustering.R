library(cluster)

#' kmedoids
#' 
cluster.kmedoids <- function(distance_mat, n_clusters) {
  output <- vector("list", length(n_clusters))
  for (i in seq_along(n_clusters)) {
    output[[i]] <- pam(distance_mat, k = n_clusters[i], diss=TRUE,
                       nstart = 10,
                       cluster.only = TRUE)
  }
  output
}
