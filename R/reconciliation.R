
toMatrix <- function(x) {
  if (is.null(dim(x))) {
    return(matrix(x, nrow = 1))
  }
  x
}

reconcile.all <- function(S, basef, resid, methods = c("ols", "wlss", "mint", "wlsv")){
  # stopifnot(!is.null(x$basef))
  output <- list()
  idx <- c(which(rowSums(S) > 1), (NROW(S) - NCOL(S) + 1):NROW(S))
  C <- S[which(rowSums(S) > 1),,drop=FALSE]
  basef <- basef[,idx]
  resid <- resid[,idx]
  n <- length(idx)
  m <- NCOL(S)
  
  for (method in methods) {
    reconcile.method <- get(paste0("reconcile.", method))
    output[[method]] <- unname(toMatrix(reconcile.method(C, basef, resid)[, c(1, (n-m+1):n)]))
  }
  output
}

reconcile.ols <- function(C, basef, resid) {
  # rm <- S %*% solve(t(S) %*% S, t(S) )
  # basef %*% t(rm)
  FoReco::htsrec(basef, "ols", C=C)$recf
}

reconcile.mint <- function(C, basef, resid) {
  FoReco::htsrec(basef, "shr", C=C, res = resid)$recf
}

reconcile.wlsv <- function(C, basef, resid) {
  FoReco::htsrec(basef, "wls", C=C, res = resid)$recf
}

reconcile.wlss <- function(C, basef, resid) {
  FoReco::htsrec(basef, "struc", C=C, res = resid)$recf
}

lambda_estimate <- function(error){
  timeT = dim(error)[1]
  covm = t(error) %*% error/timeT
  xs = apply(error, 1, function(x){x / sqrt(diag(covm))}) %>% t()
  corm = t(xs) %*% xs / timeT
  diag(corm) = 0
  d = sum(corm^2)
  xs2 = xs^2
  v = 1/(timeT*(timeT-1))*(t(xs2) %*% xs2 - 1/timeT*(t(xs) %*% xs)^2)
  diag(v) = 0
  lamb <- max(min(c(sum(v)/d, 1)), 0)
  
  lamb * diag(diag(covm)) + (1-lamb) * covm
}
