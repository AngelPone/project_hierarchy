


reconcile.wlss <- function(x){
  stopifnot(!is.null(x$basef))
  
  S <- rbind(x$S, do.call(rbind, lapply(x$nl, function(g) { g$S })))
  basef <- cbind(x$basef, do.call(cbind, lapply(x$nl, function(g){ g$basef })))
  resid <- residuals(x)
  
  W <- diag(1/apply(resid, 2, var))
  
  rm <- S %*% solve(t(S) %*% W %*% S, t(S) %*% W)
  
  rf <- basef %*% t(rm)
  x$rf <- rf[,1:NROW(x$S)]
  x
}
