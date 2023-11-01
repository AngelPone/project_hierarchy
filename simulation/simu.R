library(foreach)
library(forecast)
library(dplyr)
library(ggplot2)
source("R/reconciliation.R")

set.seed(1031)

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

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






simulate.forecast <- function(series) {

  
  bts <- t(series[,1:(NCOL(series) - 12)])
  tts <- t(series[,(NCOL(series) - 11):NCOL(series)])
  
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
  nl[[1]] <- grp2S(grp)
  for (i in 2:51) {
    nl[[i]] <- grp2S(sample(grp, length(grp)))
  }
  
  all_S <- rbind(rep(1,80), do.call(rbind, nl), diag(80))
  allts <- bts %*% t(all_S)
  bf <- foreach(x=iterators::iter(allts, by="column"), .packages = "forecast") %dopar% {
    mdl <- ets(ts(x, frequency = 12))
    fcasts <- forecast(mdl, h=12)
    list(fcasts = as.numeric(fcasts$mean), resid = as.numeric(residuals(mdl, type = "response")))
  }
  
  recfs <- list()
  for (i in seq_along(nl)) {
    C <- rbind(rep(1, 80), nl[[i]])
    basef <- do.call(cbind, lapply(bf[c(1, (4*i-2):(4*i+1), 206:285)], function(x) {x$fcasts}))
    resid <- do.call(cbind, lapply(bf[c(1, (4*i-2):(4*i+1), 206:285)], function(x) {x$resid}))
    recfs[[i]] <- reconcile.mint(C, basef, resid)[,c(1, 6:85)]
  }
  
  C <- matrix(1, ncol=80)
  basef <- do.call(cbind, lapply(bf[c(1, 206:285)], function(x) {x$fcasts}))
  resid <- do.call(cbind, lapply(bf[c(1, 206:285)], function(x) {x$resid}))
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
  methods <- c("no-cluster", "correct-cluster", "random-10", "random-20", "random-50")
  
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

FREQUENCIES <- c(12, 4, 12, 4)

for (i in 1:121) {
  print(sprintf("%s %s", Sys.time(), i))
  simulated_series <- rbind(
    simulate.season(freq = FREQUENCIES[1], length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = FALSE),
    simulate.season(freq = FREQUENCIES[2], length = 12, nseries = 20, allow.trend = TRUE, oppsite.season = TRUE),
    simulate.season(freq = FREQUENCIES[3], length = 12, nseries = 20, allow.trend = FALSE, oppsite.season = FALSE),
    simulate.season(freq = FREQUENCIES[4], length = 12, nseries = 20, allow.trend = FALSE, oppsite.season = TRUE)
  )
  
  output <- as_tibble(simulate.forecast(simulated_series)) %>% 
    mutate(batch = i) %>%
    rbind(output)
  generated_series[[i]] <- simulated_series
}

saveRDS(list(acc = output, series = generated_series), "simulation/simulation.rds")








