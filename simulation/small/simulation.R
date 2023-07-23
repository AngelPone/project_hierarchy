N <- 1000
TIME_WINDOW <- 132
WARM_UP <- 120
TRAIN_WINDOW <- 120
FREQ <- 12
GROUPS <- c(2, 2)
ALPHA <- list(upper = c(0.95, 0.95, 0.85, 0.85), lower = c(0.9, 0.9, 0.8, 0.8))
BETA <- list(upper = c(1e-5, 1e-5, 0, 0), lower = c(1e-6, 1e-6, 0, 0))
GAMMA <- list(upper = c(0.01, 0.01, 0.005, 0.005), lower = c(0.005, 0.005, 0, 0))
SPIKE <- c(2, 2, 8, 8)
SIGMA <- diag(c(4, 4, 20, 20))
SIGMA[2, 1] <- 2
SIGMA[1, 2] <- 2
SIGMA[3, 4] <- 10
SIGMA[4, 3] <- 10
SIGMA <- SIGMA / 100

sim_param <- function(x){
  sapply(1:length(x$upper), function(y){runif(1, x$lower[y], x$upper[y])})
}

sim <- function(){
  Fmat <- rbind(c(1, 1, rep(0, FREQ)), 
                c(0, 1, rep(0, FREQ)), 
                c(rep(0, FREQ+1), 1),
                cbind(0, 0, diag(FREQ - 1), 0))
  wvec <- c(1, 1, rep(0, FREQ - 1), 1)
  alpha <- sim_param(ALPHA)
  beta <- sim_param(BETA)
  gamma <- sim_param(GAMMA)
  
  epsilon <- MASS::mvrnorm(TIME_WINDOW+TIME_WINDOW, mu=rep(0, sum(GROUPS)), Sigma = SIGMA)
  
  state <- lapply(1:sum(GROUPS), function(x){
    l0 <- rnorm(1)
    b0 <- runif(1, 0.00001, 0.00002)
    s0 <- append(rnorm(FREQ - 1, sd=0.5), rnorm(1, 2, 0.01), SPIKE[x] - 1)
    s0 <- s0 - mean(s0)
    s0 <- s0[FREQ : 1]
    state <- matrix(c(l0, b0, s0), nrow = 1)
    gvec <- c(alpha[x], beta[x], gamma[x], rep(0, FREQ - 1))
    for (t in 1:(TIME_WINDOW + WARM_UP)){
      state <- rbind(state, as.vector(Fmat %*% state[t,]) + gvec * epsilon[t, x])
    }
    state[(WARM_UP + 1):(WARM_UP + TIME_WINDOW),]
  })
  y <- sapply(seq_along(state), function(x){
    as.vector(state[[x]] %*% wvec + epsilon[(WARM_UP+1): (WARM_UP + TIME_WINDOW), x])
  }, simplify = "array")
  list(y = y, params=list(alpha=alpha, beta=beta, gamma=gamma, epsilon=epsilon, state=state))
}



plot_series <- function(x){
  x <- x$y
  par(mfrow=c(3, 3))
  plot(rowSums(x), type='l', main="total", ylab = "y")
  plot(rowSums(x[,1:2]), type="l", main="middle 1", ylab = "y")
  plot(rowSums(x[,3:4]), type="l", main="middle 2", ylab = "y")
  for (i in 1:NCOL(x)){
    plot(x[,i], type='l', main=paste0("Series ", i), ylab = "y")
  }
}

produce_basef <- function(obj){
  f <- function(x){
    mdl <- forecast(ets(ts(x, frequency = 12), model="AAA", damped=FALSE, additive=TRUE), h=12)
    c(as.numeric(fitted(mdl)), mdl$mean)
  }
  obj$f <- list()
  obj$f$bottom <- sapply(1:4, function(x){ f(obj$y[1:TRAIN_WINDOW, x]) },
                         simplify = "array")
  obj$f$total <- f(rowSums(obj$y)[1:TRAIN_WINDOW])
  obj
}

library(dplyr)
library(foreach)
set.seed(42)

objs <- lapply(1:N, function(x){sim()})


cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

objs <- foreach(d=iterators::iter(objs), .errorhandling="pass", 
                .packages=c("forecast")) %dopar%
  {produce_basef(d)}

saveRDS(objs, "data/objs.rds")

parallel::stopCluster(cl)


