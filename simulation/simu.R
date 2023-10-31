library(foreach)
library(forecast)
library(dplyr)
source("R/reconcilation.R")

simulate.season <- function(freq = 12, length = 12, nseries = 50, oppsite.season = TRUE, allow.trend = TRUE, damped = FALSE) {
  seasons <- do.call(rbind,
          lapply(1:nseries, function(x) {
            a <- c(runif(freq / 2, 0, 1), runif(freq / 2, 2, 3))
            od <- rep(1:(freq/2), each = 2)
            if (oppsite.season) {
              od[seq(2, freq, 2)] <- od[seq(2, freq, 2)] + freq/2
            } else {
              od[seq(1, freq, 2)] <- od[seq(1, freq, 2)] + freq/2
            }
            rep(a[od], length * length / freq)
          }))
  
  trend <- 0
  if (allow.trend) {
    trend <- seq(0, by = 0.001, length = length * freq * nseries * length / freq) +
      rnorm(length * freq * nseries * length / freq, 0.005)
  }
  
  errors <- matrix(rnorm(length * freq * nseries * length / freq, sd=0.5), nseries)
  trend + errors + seasons
}


FREQUENCIES <- c(12, 4, 12, 4)
visualizeSimu <- function(series) {
  grpidx <- c("t-s-12", "t-os-4", "s-12", "os-4")
  grps <- rep(grpidx, each=20)
  pca <- prcomp(simulated_series, scale.=TRUE)
  
  par(mfrow = c(2,2))
  for (i in 0:3) {
    plot(ts(series[20*i+1,], frequency = FREQUENCIES[i+1]), ylab = "series", main = grpidx[i+1])
  }
  ggplot(mapping = aes(x = pca$x[,1], y=pca$x[,2], color = grps, group=grps)) +
    geom_point()
  
}


visualizeSimu(simulated_series)

simulate.forecast <- function(series) {

  
  bts <- t(series[,1:(NCOL(series) - 12)])
  tts <- t(series[,(NCOL(series) - 12):NCOL(series)])
  
  grp2S <- function(grp) {
    do.call(rbind, lapply(unique(grp), function(x) {
      S_row <- vector("integer", length(grp))
      S_row[which(grp == x)] <- 1
      S_row
    }))
  }
  # cluster
  nl <- vector("list", 51)
  grp <- rep(1:4, each = 20)
  nl[[1]] <- list(S = grp2S(grp))
  for (i in 2:51) {
    nl[[i]] <- list(S = grp2S(sample(grp, length(grp))))
  }
  
  all_S <- rbind(rep(1,80), do.call(rbind, nl), diag(80))
  allts <- bts %*% t(all_S)
  bf <- foreach(x=iterators::iter(allts, by="column")) %dopar% {
    mdl <- ets(ts(x, frequency = 12))
    fcasts <- forecast(mdl, h=12)
    list(fcasts = as.numeric(fcasts$mean), resid = residuals(mdl, type = "response"))
  }
  
  for (i in seq_along(nl)) {
    S <- nl[[i]]
    C <- rbind(rep(1, 80), S)
    basef <- do.call(cbind, lapply(bf[c(1, (4*i-2):(4*i+2)), 206:285], function(x) {x$fcasts}))
    resid <- do.call(cbind, lapply(bf[c(1, (4*i-2):(4*i+2)), 206:285], function(x) {x$resid}))
    nl[[i]]$recf <- reconcile.mint(C, basef, resid)[,c(1, 6:85)]
  }
  
  
  C <- matrix(1, ncol=80)
  basef <- do.call(cbind, lapply(bf[c(1, 206:285)], function(x) {x$fcasts}))
  resid <- do.call(cbind, lapply(bf[c(1, 206:285)], function(x) {x$resid}))
  recfs <- list()
  recf[[1]] <- reconcile.mint(C, basef, resid) # no cluster
  recf[[2]] <- nl[[1]]$recf # correct cluster
  # random 10, 20, 50
  recf[[3]] <- apply(do.call(abind::abind, lapply(nl[2:51], function(x) x$recf))[1:10,,],
                     c(2, 3), mean)
  recf[[4]] <- apply(do.call(abind::abind, lapply(nl[2:51], function(x) x$recf))[1:20,,],
                     c(2, 3), mean)
  recf[[5]] <- apply(do.call(abind::abind, lapply(nl[2:51], function(x) x$recf))[1:50,,],
                     c(2, 3), mean)
  methods <- c("no-cluster", "correct-cluster", "random-10", "random-20", "random-50")
  rmse <- vector("list", 5)
  for (i in seq_along(recf)) {
    rmse[[i]] <- sqrt(colMeans((recf[[i]] - tts)^2))
  }
  
  list(methods=methods, rmse=rmse)
}

output <- NULL

for (i in 1:121) {
  simulated_series <- rbind(
    simulate.season(freq = FREQUENCIES[1], length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = FALSE),
    simulate.season(freq = FREQUENCIES[2], length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = TRUE),
    simulate.season(freq = FREQUENCIES[3], length = 12, nseries = 20, allow.trend = FALSE, oppsite.season = FALSE),
    simulate.season(freq = FREQUENCIES[4], length = 12, nseries = 20, allow.trend = FALSE, oppsite.season = TRUE)
  )
  
  output <- rbind(as_tibble(simulate.forecast(simulated_series)) %>% mutate(batch = i))
}

saveRDS(output, "simulate/simulation.rds")








