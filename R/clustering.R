library(cluster)

cluster.random <- function(distance_mat, n_clusters) {
  
  n <- NROW(distance_mat)
  S <- vector("list", length(n_clusters))
  for (i in seq_along(n_clusters)) {
    n_cluster <- n_clusters[i]
    groups <- rep(1:n_cluster, length = n)[sample(1:n, n)]
    S[[i]] <- do.call(rbind, lapply(1:n_cluster, function(x){
      S_row <- vector("numeric", n)
      S_row[which(groups == x)] <- 1
      S_row
    }))
  }
  S
}


#' kmedoids
#' 
cluster.kmedoids <- function(distance_mat, n_clusters) {
  output <- vector("list", length(n_clusters))
  for (i in seq_along(n_clusters)) {
    output[[i]] <- pam(distance_mat, k = n_clusters[i], diss=TRUE,
                       nstart = 10,
                       cluster.only = TRUE)
  }
  grplst2Slst <- function(grplsts) {
    lapply(grplsts, function(grpvec) {
      do.call(rbind, lapply(unique(grpvec), function(grp){
        S_row <- vector("numeric", NCOL(distance_mat))
        S_row[which(grpvec == grp)] <- 1
        S_row
      }))
    })
  }
  
  grplst2Slst(output)
}

#' hierarchical clustering
#' 
#' @param method linkage method: see ?cluster::agnes
cluster.hcluster <- function(distance_mat, method) {
  hc <- agnes(distance_mat, diss = FALSE, method = method,
        keep.data = FALSE, keep.diss = FALSE)
  
  S <- matrix(0, NROW(hc$merge) - 1, NROW(distance_mat))
  
  for (i in 1:NROW(S)) {
    cur_idx <- hc$merge[i,]
    for (idx in cur_idx) {
      if (idx < 0) {
        S[i, abs(idx)] <- 1
      } else {
        S[i, which(S[idx, ] == 1)] <- 1
      }
    }
  }
  list(S)
}



