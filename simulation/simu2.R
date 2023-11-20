library(foreach)
library(forecast)
library(dplyr)
library(ggplot2)
source("R/reconciliation.R")

set.seed(43)

cl <- parallel::makeCluster(3)
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
  for (i in 2:51) {
    nl[[i]] <- grp2S(sample(grp, length(grp)))
  }
  nl[[52]] <- grp2S(rep(1:3, each=40))
  nl[[53]] <- grp2S(rep(c(1,2, 1, 2, 1, 2), each=20))
  nl[[54]] <- grp2S(rep(c(1,2, 1), each=40))
  
  all_S <- rbind(rep(1,m), do.call(rbind, nl))
  allts <- bts %*% t(all_S)
  bf <- foreach(x=iterators::iter(allts, by="column"), .packages = "forecast") %dopar% {
    mdl <- ets(ts(x, frequency = 12))
    fcasts <- forecast(mdl, h=12)
    list(fcasts = as.numeric(fcasts$mean), resid = as.numeric(residuals(mdl, type = "response")))
  }
  bf_bottom <- lapply(1:m, function(x){
    fcasts <- mean(bts[,x])
    list(fcasts = rep(fcasts, 12), resid = fcasts - bts[,x])
  })
  bf <- c(bf, bf_bottom)
  
  bottom_idx <- (NROW(all_S)+1):(NROW(all_S)+m)
  
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
  recf_output[[2]] <- recfs[[1]] # correct cluster
  # random 10, 20, 50
  recf_output[[3]] <- apply(do.call(abind::abind, list(recfs[2:51], along=0))[1:10,,],
                            c(2, 3), mean)
  recf_output[[4]] <- apply(do.call(abind::abind, list(recfs[2:51], along=0))[1:20,,],
                            c(2, 3), mean)
  recf_output[[5]] <- apply(do.call(abind::abind, list(recfs[2:51], along=0))[1:50,,],
                            c(2, 3), mean)
  recf_output[[6]] <- recfs[[52]]
  recf_output[[7]] <- recfs[[53]]
  recf_output[[8]] <- recfs[[54]]
  recf_output[[9]] <- apply(do.call(abind::abind, list(recfs[c(1, 52, 53, 54)], along=0)),
                            c(2, 3), mean)
  recf_output[[10]] <- do.call(cbind, lapply(bf[c(1, bottom_idx)], function(x) {x$fcasts}))
  methods <- c("no-cluster", "correct-cluster", "random-10", "random-20", "random-50",
               "trend-cluster", "season-cluster", "trend-cluster2", "cluster-average", "base")
  
  # evaluate
  rmse <- vector("list", 5)
  tts <- cbind(rowSums(tts), tts)
  for (i in seq_along(recf_output)) {
    rmse[[i]] <- unname(sqrt(colMeans((recf_output[[i]] - tts)^2)))
  }
  list(methods=methods, rmse=rmse)
}

output <- NULL
generated_series <- vector("list", 121)


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
  
  output <- as_tibble(simulate.forecast(simulated_series)) %>% 
    mutate(batch = i) %>%
    rbind(output)
  generated_series[[i]] <- simulated_series
}

saveRDS(list(acc = output, series = generated_series), "simulation/simulation2.rds")








