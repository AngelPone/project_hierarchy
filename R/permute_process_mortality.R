
nl2tibble <- function(x) {
  output <- vector("list", 6)
  names(output) <- c("representor", "cluster", "distance", "S", "rf", "other")
  
  for (i in seq_along(x)) {
    for (n in names(x[[i]])) {
      if (is.character(x[[i]][[n]])) {
        output[[n]] <- append(output[[n]], x[[i]][[n]])
      } else {
        output[[n]] <- append(output[[n]], list(x[[i]][[n]]))
      }
    }
  }
  tibble::as_tibble(output)
}

for (i in 0:119) {
  print(i)
  dt <- readRDS(sprintf("tourism/ets/batch_%s.rds", i))
  cluster <- sapply(dt$nl, function(x){x$cluster})
  idx_to_remove <- which(startsWith(cluster, "permute-hcluster"))[101:600]
  dt$nl <- dt$nl[-idx_to_remove]
  stopifnot(length(dt$nl) == 964)
  saveRDS(dt, sprintf("tourism/ets/batch_%s.rds", i))
}
