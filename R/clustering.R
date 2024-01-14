library(cluster, quietly = TRUE)

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

cluster.nestedkmedoids <- function(distance_mat, n_clusters) {
  
  output <- list()
  output[[1]] <- list(1:NROW(distance_mat))
  
  for (i in 1+seq_along(n_clusters)) {
    output[[i]] <- list()
    k <- 1
    for (j in seq_along(output[[i-1]])) {
      if (length(output[[i-1]][[j]]) <= n_clusters[i-1]) {
        next
      }
      res <- pam(distance_mat[output[[i-1]][[j]], output[[i-1]][[j]]], k = n_clusters[i-1], diss=TRUE,
          nstart = 10,
          cluster.only = TRUE)
      for (ij in 1:n_clusters[i-1]) {
        output[[i]][[k]] <- output[[i-1]][[j]][which(res == ij)]
        k <- k+1
      }
    }
  }
  
  list(do.call(rbind, lapply(output[2:length(output)], function(x) {
    do.call(rbind, lapply(x, function(g) {
      S_row <- vector("numeric", NROW(distance_mat))
      S_row[g] <- 1
      S_row
    }))
  })))
}


#' kmedoids
#' 
cluster.kmedoids <- function(distance_mat, n_clusters) {
  
  # calculate silhouette and determine optimal number of clusters
  max.avgwidths <- 0
  max.n_cluster <- 1
  max.pr <- NULL
  for (n_cluster in n_clusters) {
    if (n_cluster == 1) next
    pr <- pam(distance_mat, k=n_cluster, diss=TRUE)
    clus_silwidth <- summary(silhouette(pr))
    # if (min(clus_silwidth$clus.avg.widths) < 0.1) next
    if (clus_silwidth$avg.width > max.avgwidths) {
      max.avgwidths <- clus_silwidth$avg.width
      max.n_cluster <- n_cluster
      max.pr <- pr
    }
  }
  if (max.n_cluster == 1) { return(NULL) }
  
  grp2S <- function(grp) {
    grpvec <- grp$clustering
    # bad_groups <- which(summary(silhouette(grp))$clus.avg.width < 0.05)
    do.call(rbind, lapply(unique(grpvec), function(grp){
      # if (grp %in% bad_groups) {
      #   return(NULL)
      # }
      S_row <- vector("numeric", NCOL(distance_mat))
      S_row[which(grpvec == grp)] <- 1
      S_row
    }))
  }

  list(S=grp2S(max.pr), 
       info = list(n_cluster = max.n_cluster, 
                   width = summary(silhouette(max.pr))$clus.avg.widths,
                   avgwidth = summary(silhouette(max.pr))$avg.width,
                   medoids = max.pr$id.med)
       )
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



