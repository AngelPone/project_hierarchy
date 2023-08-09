representator.ts <- function(x) {
  apply(x$bts, 2, function(x){
    x / sd(x)
  })
}


representator.error <- function(x){
  stopifnot(!is.null(x$basef))
  n <- NROW(x$S)
  m <- NCOL(x$S)
  apply(x$resid[,(n-m+1):n], 2, function(x) {
    x / sd(x)
  })
}
