

S <- rbind(c(1,1,1), c(0, 1, 1), diag(3))


weights <- function(W) {
  W_inv <- solve(W)
  S %*% solve(t(S) %*% W_inv %*% S, t(S) %*% W_inv)
}


weights_plot <- function(offdiag, title = "Setting 1") {
  
  
  construct_W <- function(offdiag, x) {
    m <- matrix(0, 5, 5)
    m[lower.tri(m)] <- offdiag
    m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
    diag(m) <- c(1, 1, 1, 1, 1)
    m[4, 5] <- x
    m[5, 4] <- x
    diag(sqrt(c(3, 2.1, 0.9, 1.05, 1.04))) %*% m %*% diag(sqrt(c(3, 2.1, 0.9, 1.05, 1.04)))
  }
  print(construct_W(offdiag, 0.5))
  
  weights_all <- vector("list", 5)
  
  for (i in seq(-0.99, 0.99, by = 0.01)) {
    W <- construct_W(offdiag, i)
    weight_i <- weights(W)
    for (j in 1:5){
      weights_all[[j]] <- rbind(weights_all[[j]], c(i, weight_i[,j]))
    }
  }
  library(dplyr)
  for (j in 1:5) {
    colnames(weights_all[[j]]) <- c("cor", "A", "B", "C", "D", "E")
    weights_all[[j]] <- as_tibble(weights_all[[j]]) %>%
      mutate(weights_col = c("A", "B", "C", "D", "E")[j])
  }
  
  weights_all <- do.call(rbind, weights_all)
  
  weights_all <- weights_all %>% tidyr::pivot_longer(names_to = "weights_row", cols = c("A", "B", "C", "D", "E"))
  
  library(ggplot2)
  
  ggplot(weights_all)  +
    geom_line(mapping = aes(x = cor, y = value, color = weights_row, group=weights_row)) +
    facet_wrap(~ weights_col) +
    ggtitle(title)
}


# # Setting 1
# weights_plot(0)
# 
# # Setting 2
# 
# offdiag <- runif(10, 0, 0.5)
# weights_plot(offdiag, title = "Setting 2")
# 
# # Setting 3
# offdiag <- runif(10, -0.5, 0)
# weights_plot(offdiag, title = "Setting 3")
# 
# # Setting 4
# offdiag <- runif(10, -0.3, 0.3)
# weights_plot(offdiag, title = "Setting 4")






# AR(1) errors
library(forecast)
library(dplyr)
source("R/reconciliation.R")
step2 <- function(epsilon_cor, epsilon_var, ma_theta) {
  Wb <- diag(3)
  Wb[lower.tri(Wb)] <- c(0, 0, epsilon_cor)
  Wb[upper.tri(Wb)] <- t(Wb)[upper.tri(t(Wb))]
  
  Wb <- diag(sqrt(epsilon_var)) %*% Wb %*% diag(sqrt(epsilon_var))
  errors <- MASS::mvrnorm(1001, rep(0, 3), Wb)
  AR <- diag(c(0, ma_theta, ma_theta))
  # AR[2,3] <- offdiag[1]
  # AR[3,2] <- offdiag[2]
  series <- matrix(0, 1001, 3)
  for (i in 2:1001) {
    series[i,] <- AR %*% errors[i-1, ] + errors[i,]
  }
  series <- series[602:1001,]
  
  # print(cov(series))
  series[,1] <- seq(from=1, to=40, length = 400) + series[,1]
  series[,2] <- seq(from=1, to=36, length = 400) + series[,2]
  series[,3] <- seq(from=1, to=-36, length = 400) + series[,3]
  
  # plot(ts(series))
  
  S2 <- S
  S2[2,] <- c(1, 1, 0)
  
  series2 <- series %*% t(S2)
  series <- series %*% t(S)
  # Cross Validation
  init_window <- 300
  step <- 10
  
  mat2list <- function(x){
    r <- lapply(seq_len(ncol(x)), function(i) unname(x[,i]))
    names(r) <- colnames(x)
    r
  } 
  accs <- vector("list", 4)
  names(accs) <- c("base", "error", "ts", "noclustering")
  i <- 10
  train <- series[1:(290+10*i),]
  test <- series[(290+10*i+1):(300+10*i),]
  test2 <- series2[(290+10*i+1):(300+10*i),]
    
  fcasts <- list()
  for (j in 1:5) {
    if (j <= 2) {
      mdl <- forecast(auto.arima(train[, j], seasonal = FALSE), h=10)
    } else {
      mdl <- forecast(arima(train[, j], order = c(0, 1, 0)), h=10)
    }
    fcasts[[j]] <- list(resid = as.numeric(residuals(mdl)), fcasts = as.numeric(mdl$mean))
  }
  
  resids <- do.call(cbind, lapply(fcasts, function(x){
    x$resid
  }))
  fcasts <- do.call(cbind, lapply(fcasts, function(x){
    x$fcasts
  }))
  
  fcasts_B2 <- forecast(auto.arima(series2[1:(290+10*i), 2]), h=10)
  resids2 <- resids
  resids2[,2] <- as.numeric(residuals(fcasts_B2))
  fcasts2 <- fcasts
  fcasts2[,2] <- as.numeric(fcasts_B2$mean)
  
  reconf3 <- reconcile.mint(S[c(1, 3, 4, 5),], fcasts[,c(1, 3, 4, 5)], resids[,c(1, 3, 4, 5)])
  reconf <- reconcile.mint(S, fcasts, resids)
  reconf2 <- reconcile.mint(S2, fcasts2, resids2)
  
  rmse <- function(x, y) { 
    r <- sapply(c(1, 2, 3, 4, 5), function(i){
      sqrt(mean((x[,i]-y[,i])^2))
    })
    r
  }
  
  rmse_noclustering <- function(x, y) {
    r <- sapply(c(1, 2, 3, 4), function(i){
      sqrt(mean((x[,i]-y[,i])^2))
    })
    r <- c(r[1], NA, r[2:4])
    r
  }
  output <- data.frame(base = rmse(fcasts, test), error = rmse(reconf, test),
    ts = rmse(reconf2, test2), noclustering = rmse_noclustering(reconf3, test2[,c(1,3,4,5)]))
  output$name <- c("A", "B", "C", "D", "E")
  output
}

r1 <- NULL
for (i in 1:100) {
  r1 <- step2(epsilon_cor = 0.9,epsilon_var = c(100, 4, 4), ma_theta = 0.3) %>%
    mutate(i = .env$i) %>%
    rbind(r1)
}

r2 <- NULL
for (i in 1:100) {
  r2 <- step2(epsilon_cor = 0,epsilon_var = c(100, 4, 4), ma_theta = 0.3) %>%
    mutate(i = .env$i) %>%
    rbind(r2)
}

r3 <- NULL
for (i in 1:100) {
  r3 <- step2(epsilon_cor = -0.9, epsilon_var = c(100, 4, 4), ma_theta = 0.3) %>%
    mutate(i = .env$i) %>%
    rbind(r3)
}
# 
# 
r1 %>% mutate(error = error/base, ts = ts/base, noclustering=noclustering/base) %>% group_by(name) %>% summarise_at(c("error", "ts", "noclustering"), mean)
r2 %>% mutate(error = error/base, ts = ts/base, noclustering=noclustering/base) %>% group_by(name) %>% summarise_at(c("error", "ts", "noclustering"), mean)
r3 %>% mutate(error = error/base, ts = ts/base, noclustering=noclustering/base) %>% group_by(name) %>% summarise_at(c("error", "ts", "noclustering"), mean)

r3 %>% group_by(name) %>% summarise_at(c("error", "ts", "noclustering", "base"), mean)
