

compute_entropy <- function(dataset, bf_method) {
  files <- dir(paste0(dataset, "/", bf_method))
  batches <- files[which(startsWith(files, "store_"))]
  
  output <- list()
  for (batch in batches) {
    dt <- readRDS(paste0(dataset, "/", bf_method, "/", batch))
    batch_n <- as.integer(strsplit(strsplit(batch, "_")[[1]][[2]], ".", fixed = TRUE)[[1]][[1]])
    error <- dt$data$resid
    output[[batch_n+1]] <- apply(error, 2, tsfeatures::entropy)
  }
  do.call(rbind, output)
}

compute_tsentropy <- function(dataset, bf_method) {
  files <- dir(paste0(dataset, "/", bf_method))
  batches <- files[which(startsWith(files, "store_"))]
  
  output <- list()
  for (batch in batches) {
    dt <- readRDS(paste0(dataset, "/", bf_method, "/", batch))
    batch_n <- as.integer(strsplit(strsplit(batch, "_")[[1]][[2]], ".", fixed = TRUE)[[1]][[1]])
    bts <- dt$data$bts
    output[[batch_n+1]] <- apply(bts, 2, function(x){
      tsfeatures::entropy(ts(x, frequency = 12))
    })
  }
  do.call(rbind, output)
}

mortality_ets <- compute_entropy("mortality", "ets")
mortality_arima <- compute_entropy("mortality", "arima")


tourism_ets <- compute_entropy("tourism", "ets")
tourism_arima <- compute_entropy("tourism", "arima")
