TRAIN_WINDOW <- 120
TIME_WINDOW <- 132

S1 <- rbind(rep(1, 4), c(1, 1, 0, 0), c(0, 0, 1, 1), diag(4))
S2 <- rbind(rep(1, 4), c(1, 0, 1, 0), c(0, 1, 0, 1), diag(4))
S3 <- rbind(rep(1, 4), c(1, 0, 0, 1), c(0, 1, 1, 0), diag(4))

produce_mid_basef <- function(obj){
  f <- function(x){
      mdl <- forecast(ets(ts(x, frequency = 12), model="AAA", damped=FALSE, additive=TRUE), h = 12)
      c(fitted(mdl), mdl$mean)
    }
  obj$f$middle <- list(apply(obj$y[1:TRAIN_WINDOW,] %*% t(S1[2:3,]), 2, f),
                       apply(obj$y[1:TRAIN_WINDOW,] %*% t(S2[2:3,]), 2, f),
                       apply(obj$y[1:TRAIN_WINDOW,] %*% t(S3[2:3,]), 2, f))
  obj
}

lambda_estimate <- function(error){
  timeT <- dim(error)[1]
  covm <- cov(error)
  xs <- t(apply(error, 1, function(x){x / sqrt(diag(covm))}))
  corm <- stats::cov2cor(covm)
  diag(corm) <- 0
  d <- sum(corm^2)
  xs2 <- xs^2
  v <- 1/(timeT*(timeT-1))*(t(xs2) %*% xs2 - 1/timeT*(t(xs) %*% xs)^2)
  diag(v) <- 0
  max(min(c(sum(v)/d, 1)), 0)
}

cal_recMat <- function(S, residuals){
  # wlss
  wlss_weights <- diag(1/rowSums(S))
  # wlsv
  cov_mat <- cov(residuals)
  wlsv_weights <- diag(1/diag(cov_mat))
  # shrinkage
  lamb <-lambda_estimate(residuals)
  shrink_weights <- (1-lamb) * cov_mat + lamb * wlsv_weights
  shrink_weights <- solve(shrink_weights)
  
  return(list(ols = solve(t(S) %*% S, t(S)),
              wlss = solve(t(S) %*% wlss_weights %*% S, t(S) %*% wlss_weights),
              wlsv = solve(t(S) %*% wlsv_weights %*% S, t(S) %*% wlsv_weights),
              shrinkage = solve(t(S) %*% shrink_weights %*% S, t(S) %*% shrink_weights)))
}


reconcile <- function(obj){
  S <- list(S1, S2, S3)
  obj$rf <- lapply(1:3, function(x){
    f <- obj$f
    basef <- cbind(f$total, f$middle[[x]], f$bottom)
    y <- obj$y %*% t(S[[x]])
    resid <- (basef - y)[1:TRAIN_WINDOW, ]
    recMats <- cal_recMat(S[[x]], resid)
    lapply(recMats, function(g){
      (basef[(TRAIN_WINDOW+1):TIME_WINDOW,] %*% t(S[[x]] %*% g))[,c(1, 4:7)]
    })
  })
  obj
}



library(foreach)

objs <- readRDS("data/objs.rds")

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

objs <- foreach(d=iterators::iter(objs), .errorhandling="pass", 
                .packages=c("forecast")) %dopar%
  {reconcile(produce_mid_basef(d))}

saveRDS(objs, "data/objs.rds")

parallel::stopCluster(cl)


