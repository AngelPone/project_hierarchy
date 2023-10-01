

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


# Setting 1
weights_plot(0)

# Setting 2

offdiag <- runif(10, 0, 0.5)
weights_plot(offdiag, title = "Setting 2")

# Setting 3
offdiag <- runif(10, -0.5, 0)
weights_plot(offdiag, title = "Setting 3")

# Setting 4
offdiag <- runif(10, -0.3, 0.3)
weights_plot(offdiag, title = "Setting 4")






# AR(1) errors
library(forecast)
source("R/reconciliation.R")
step2 <- function(offdiag) {
  
  Wb <- diag(3)
  Wb[lower.tri(Wb)] <- offdiag
  Wb[upper.tri(Wb)] <- t(Wb)[upper.tri(t(Wb))]
  
  Wb <- diag(sqrt(c(100, 50, 20))) %*% Wb %*% diag(sqrt(c(100, 50, 20)))
  
  errors <- MASS::mvrnorm(1000, rep(0, 3), Wb)
  
  AR <- diag(c(0, 0.3, 0.4))
  series <- matrix(0, 1001, 3)
  for (i in 2:1001) {
    series[i,] <- AR %*% series[i-1, ] + errors[i-1,]
  }
  series <- series[102:1001,]
  series[,1] <- seq(from=1, to=200, length = 900) + series[,1]
  series[,2] <- seq(from=1, to=180, length = 900) + series[,2]  
  
  plot(ts(series))
  
  S2 <- S
  S2[2,] <- c(1, 1, 0)
  
  series2 <- series %*% t(S2)
  series <- series %*% t(S)
  # Cross Validation
  init_window <- 800
  step <- 10
  
  mat2list <- function(x){
    r <- lapply(seq_len(ncol(x)), function(i) unname(x[,i]))
    names(r) <- colnames(x)
    r
  } 
  accs <- vector("list", 3)
  names(accs) <- c("base", "mint", "mint2")
  for (i in 1:10) {
    train <- series[1:(790+10*i),]
    test <- series[(791+10*i):(800+10*i),]
    
    fcasts <- lapply(iterators::iter(train, by="col"), function(s){
      mdl <- forecast(auto.arima(s, max.p = 1, max.order=2, seasonal = FALSE), h=10)
      list(resid = as.numeric(residuals(mdl)), fcasts = as.numeric(mdl$mean))
    })
    
    resids <- do.call(cbind, lapply(fcasts, function(x){
      x$resid
    }))
    fcasts <- do.call(cbind, lapply(fcasts, function(x){
      x$fcasts
    }))
    
    fcasts_B2 <- forecast(auto.arima(series2[1:(790+10*i), 2], seasonal = FALSE), h=10)
    resids2 <- resids
    resids2[,2] <- as.numeric(residuals(fcasts_B2))
    fcasts2 <- fcasts
    fcasts2[,2] <- as.numeric(fcasts_B2$mean)
    
    
    reconf <- reconcile.mint(S, fcasts, resids)
    reconf2 <- reconcile.mint(S2, fcasts2, resids2)
    
    rmse <- function(x, y) { 
      r <- sapply(c(1, 3, 4, 5), function(i){
        sqrt(mean((x[,i]-y[,i])^2))
      })
      names(r) <- c("A", "C", "D", "E")
      r
    }
    accs$base <- rbind(accs$base, rmse(fcasts, test))
    accs$mint <- rbind(accs$mint, rmse(reconf, test))
    accs$mint2 <- rbind(accs$mint2, rmse(reconf2, test))
  }
  print(cov(resids))
  accs$base <- mat2list(accs$base)
  accs$mint <- mat2list(accs$mint)
  accs$mint2 <- mat2list(accs$mint2)
  as_tibble(accs) %>% mutate(names = c("A", "C", "D", "E")) %>% tidyr::unnest(cols = c("base", "mint", "mint2"))
}



r1 <- step2(c(0, 0, 0.9))
r1 %>% mutate(mint = mint/base, mint2 = mint2/base) %>% group_by(names) %>% summarise_at(c("mint", "mint2"), mean)


r2 <- step2(c(0, 0, 0))
r2 %>% mutate(mint = mint/base, mint2 = mint2/base) %>% group_by(names) %>% summarise_at(c("mint", "mint2"), mean)

r3 <- step2(c(0, 0, -0.9))
r3 %>% group_by(names) %>% summarise_at(c("base", "mint", "mint2"), mean)
r3 %>% mutate(mint = mint/base, mint2 = mint2/base) %>% group_by(names) %>% summarise_at(c("mint", "mint2"), mean)





