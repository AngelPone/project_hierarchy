
toMatrix <- function(x) {
  if (is.null(dim(x))) {
    return(matrix(x, nrow = 1))
  }
  x
}

reconcile.all <- function(S, basef, resid, methods = c("ols", "wlss", "mint", "wlsv")){
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
