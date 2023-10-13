entropy <- function(x, order = NULL) {
  spec <- try(stats::spec.ar(na.contiguous(x), plot=FALSE, method='burg',
                             n.freq = ceiling(length(x)/2 + 1), order = order))
  if ("try-error" %in% class(spec)) {
    entropy <- NA
  } else {
    fx <- c(rev(spec$spec[-1]),spec$spec)/ length(x)
    fx <- fx/sum(fx)
    prior.fx = rep(1 / length(fx), length = length(fx))
    prior.weight = 0.001
    fx <- (1 - prior.weight) * fx + prior.weight * prior.fx
    entropy <- pmin(1, -sum(fx * log(fx, base = length(x))))
  }
  return(c(entropy = entropy))
}


compute_entropy <- function(dataset, bf_method, order = NULL) {
  files <- dir(paste0(dataset, "/", bf_method))
  batches <- files[which(startsWith(files, "store_"))]
  
  output <- list()
  for (batch in batches) {
    dt <- readRDS(paste0(dataset, "/", bf_method, "/", batch))
    batch_n <- as.integer(strsplit(strsplit(batch, "_")[[1]][[2]], ".", fixed = TRUE)[[1]][[1]])
    error <- dt$data$resid
    if (is.null(order)) {
      output[[batch_n+1]] <- apply(error, 2, entropy)
    } else {
      output[[batch_n+1]] <- do.call(cbind, lapply(order, function(o) {
        apply(error, 2, entropy, order = o)
      }))
    }
    
  }
  do.call(abind::abind, list(output, along=0),)
}

# compute_tsentropy <- function(dataset, bf_method) {
#   files <- dir(paste0(dataset, "/", bf_method))
#   batches <- files[which(startsWith(files, "store_"))]
#   
#   output <- list()
#   for (batch in batches) {
#     dt <- readRDS(paste0(dataset, "/", bf_method, "/", batch))
#     batch_n <- as.integer(strsplit(strsplit(batch, "_")[[1]][[2]], ".", fixed = TRUE)[[1]][[1]])
#     bts <- dt$data$bts
#     output[[batch_n+1]] <- apply(bts, 2, function(x){
#       entropy(ts(x, frequency = 12))
#     })
#   }
#   
# }



mortality_ets <- compute_entropy("mortality", "ets", order = 1:12)
mortality_arima <- compute_entropy("mortality", "arima", order = 1:12)

tourism_ets <- compute_entropy("tourism", "ets", order = 12)
tourism_arima <- compute_entropy("tourism", "arima", order = 12)


saveRDS(list(mortality = list(ets = mortality_ets, arima = mortality_arima),
             tourism = list(ets = tourism_ets, arima = tourism_arima)), "entropy.rds")

par(mfrow = c(2,2))


hist(mortality_ets)
hist(mortality_arima)

hist(tourism_arima)
hist(tourism_ets)


