
toMatrix <- function(x) {
  if (is.null(dim(x))) {
    return(matrix(x, nrow = 1))
  }
  x
}

reconcile.all <- function(x, methods = c("ols", "wlss", "mint", "wlsv")){
  stopifnot(!is.null(x$basef))
  
  # without clusters
  if (is.null(x$rf)) {
    S <- x$S
    resid <- x$resid
    basef <- x$basef
    x$rf <- list()
    for (method in methods) {
      reconcile.method <- get(paste0("reconcile.", method))
      x$rf[[method]] <- toMatrix(reconcile.method(S, basef, resid))
    }
  }
  # with clusters
  
  for (i in seq_along(x$nl)) {
    stopifnot(!is.null(x$nl[[i]]$basef))
    nl <- x$nl[[i]]
    if (!is.null(nl$rf)) { next }
    else { x$nl[[i]]$rf <- list() }
    
    S <- rbind(x$S, nl$S)
    basef <- cbind(x$basef, nl$basef)
    resid <- cbind(x$resid, nl$resid)
    
    for (method in methods) {
      reconcile.method <- get(paste0("reconcile.", method))
      x$nl[[i]]$rf[[method]] <- toMatrix(reconcile.method(S, basef, resid)[, 1:NROW(x$S)])
    }
  }
  x
}

reconcile.ols <- function(S, basef, resid) {
  rm <- S %*% solve(t(S) %*% S, t(S) )
  basef %*% t(rm)
}

reconcile.mint <- function(S, basef, resid) {
  lamb <- lambda_estimate(resid)
  W <- cov(resid)
  W <- solve((1-lamb) * W + lamb * diag(diag(W)))
  rm <- S %*% solve(t(S) %*% W %*% S, t(S) %*% W)
  
  basef %*% t(rm)
}

reconcile.wlsv <- function(S, basef, resid) {
  W <- diag(1/apply(resid, 2, var))
  rm <- S %*% solve(t(S) %*% W %*% S, t(S) %*% W)
  basef %*% t(rm)
}

reconcile.wlss <- function(S, basef, resid) {
  W <- diag(1/rowSums(S))
  rm <- S %*% solve(t(S) %*% W %*% S, t(S) %*% W)
  basef %*% t(rm)
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
  max(min(c(sum(v)/d, 1)), 0)
}
