library(foreach)
library(forecast)
library(dplyr)
library(ggplot2)
source("R/utils.R")

set.seed(43)

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

simulate.season <- function(freq = 12, length = 12, nseries = 50, oppsite.season = TRUE, allow.trend = TRUE, opposite.trend = FALSE) {
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
    if (opposite.trend) {
      trend <- seq(0, by = -0.002, length = length * freq * nseries * length / freq) +
        rnorm(length * freq * nseries * length / freq, 0.007) 
    } else {
      trend <- seq(0, by = 0.001, length = length * freq * nseries * length / freq) +
        rnorm(length * freq * nseries * length / freq, 0.005)
    }
  }
  errors <- matrix(rnorm(length * freq * nseries * length / freq, sd=0.5), nseries)
  trend + errors + seasons
}

new_groups <- lapply(1:100, function(x){
  sample(120)
})

simulate.forecast <- function(series) {

  bts <- t(series[,1:(NCOL(series) - 12)])
  tts <- t(series[,(NCOL(series) - 11):NCOL(series)])
  m <- NROW(series)
  grp2S <- function(grp) {
    do.call(rbind, lapply(unique(grp), function(x) {
      S_row <- vector("integer", length(grp))
      S_row[which(grp == x)] <- 1
      S_row
    }))
  }
  # cluster
  nl <- vector("list", 54)
  grp <- rep(1:(m/20), each = 20)
  nl[[1]] <- grp2S(grp)
  for (i in 2:101) {
    nl[[i]] <- nl[[1]][,new_groups[[i-1]]]
  }
  
  grp <- grp2S(rep(1:3, each=40))
  for (i in 102:201) {
    nl[[i]] <- grp[,new_groups[[i-101]]]
  }
  nl[[302]] <- grp
  
  grp <- grp2S(rep(c(1,2, 1, 2, 1, 2), each=20))
  for (i in 202:301) {
    nl[[i]] <- grp[,new_groups[[i-201]]]
  }
  nl[[303]] <- grp
  
  grp <- grp2S(rep(c(1,2, 1), each=40))
  nl[[304]] <- grp
  for (i in 305:404) {
    nl[[i]] <- grp[, new_groups[[i-304]]]
  }
  
  all_S <- rbind(rep(1,m), do.call(rbind, nl), diag(m))
  allts <- bts %*% t(all_S)
  bf <- foreach(x=iterators::iter(allts, by="column"), .packages = "forecast") %dopar% {
    mdl <- ets(ts(x, frequency = 12))
    fcasts <- forecast(mdl, h=12)
    list(fcasts = as.numeric(fcasts$mean), resid = as.numeric(residuals(mdl, type = "response")))
  }
  
  bottom_idx <- (NROW(all_S)-(m-1)):NROW(all_S)
  
  recfs <- list()
  cumNROW <- 1
  for (i in seq_along(nl)) {
    C <- rbind(rep(1, m), nl[[i]])
    basef <- do.call(cbind, lapply(bf[c(1, (cumNROW+1):(cumNROW+NROW(nl[[i]])), bottom_idx)], function(x) {x$fcasts}))
    resid <- do.call(cbind, lapply(bf[c(1, (cumNROW+1):(cumNROW+NROW(nl[[i]])), bottom_idx)], function(x) {x$resid}))
    recfs[[i]] <- reconcile.mint(C, basef, resid)[,c(1, (1+NROW(C)):(NROW(C)+m))]
    cumNROW <- cumNROW+NROW(nl[[i]])
  }
  
  
  C <- matrix(1, ncol=m)
  basef <- do.call(cbind, lapply(bf[c(1, bottom_idx)], function(x) {x$fcasts}))
  resid <- do.call(cbind, lapply(bf[c(1, bottom_idx)], function(x) {x$resid}))
  recf_output <- list()
  recf_output[[1]] <- reconcile.mint(C, basef, resid) # no cluster
  recf_output[2:(length(recfs) + 1)] <- recfs
  recf_output
}

output <- list()
generated_series <- vector("list", 500)


for (i in 1:500) {
  print(sprintf("%s %s", Sys.time(), i))
  simulated_series <- rbind(
    simulate.season(freq = 12, length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = FALSE),
    simulate.season(freq = 12, length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = TRUE),
    simulate.season(freq = 12, length = 12, nseries = 20, allow.trend = FALSE, oppsite.season = FALSE),
    simulate.season(freq = 12, length = 12, nseries = 20, allow.trend = FALSE, oppsite.season = TRUE),
    simulate.season(freq = 12, length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = FALSE, opposite.trend = TRUE),
    simulate.season(freq = 12, length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = TRUE, opposite.trend = TRUE)
  )
  
  output[[i]] <- simulate.forecast(simulated_series)
  generated_series[[i]] <- simulated_series
}

saveRDS(list(acc = output, series = generated_series), "simulation/simulation.rds")







