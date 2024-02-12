args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

# cl <- parallel::makeCluster(8)
# doParallel::registerDoParallel(cl)

dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- 1
frequency <- 12
batch_length <- time_length - 96 + 12
batch_length <- time_length - 96 - forecast_horizon

metrics <- c("rmsse")
source("R/metrics.R")


hts.eval2 <- function(df, metrics, tts, bts) {
  if (is.null(dim(tts))) {
    tts <- c(sum(tts), tts)
    tts <- matrix(tts, nrow = 1)
  } else {
    tts <- cbind(rowSums(tts), tts)[1:forecast_horizon,,drop=FALSE]
  }
  
  bts <- cbind(rowSums(bts), bts)
  
  # df <- df %>% filter(!startsWith(cluster, "random")) %>%
  #   filter(cluster != "hcluster-random")
  #   # filter(!startsWith(cluster, "permute")) %>%
  #   # filter(cluster != "")
  
  for (metric in metrics) {
    accuracy_method <- get(paste0("metric.", metric))
    
    df[[metric]] <- foreach::foreach(g=iterators::iter(df$rf)) %do%  {
      c <- g[['mint']][1:forecast_horizon,,drop=FALSE]
      sapply(1:NCOL(c), function(x) { accuracy_method(tts[,x], c[,x], bts[,x]) } )
    }
  }
  df
}

nl2tibble <- function(x) {
  output <- vector("list", 6)
  names(output) <- c("representor", "cluster", "distance", "S", "rf", "other")
  
  for (i in seq_along(x)) {
    for (n in names(x[[i]])) {
      if (is.character(x[[i]][[n]])) {
        output[[n]] <- append(output[[n]], x[[i]][[n]])
      } else {
        output[[n]] <- append(output[[n]], list(x[[i]][[n]]))
      }
    }
  }
  tibble::as_tibble(output)
}

hts.evalbase <- function(dt, metrics) {
  if (is.null(dim(dt$tts))){
    tts <- c(sum(dt$tts), dt$tts)
    tts <- matrix(tts, nrow=1)
  } else {
    tts <- cbind(rowSums(dt$tts), dt$tts)[1:forecast_horizon,,drop=FALSE]
  }
  
  bts <- cbind(rowSums(dt$bts), dt$bts)
  output <- list()
  for (metric in metrics) {
    accuracy_method <- get(paste0("metric.", metric))
    output[[metric]] <-
      list(sapply(1:NCOL(bts), function(x) { accuracy_method(tts[,x], dt$basef[1:forecast_horizon,x], bts[,x]) } ))
  }
  output
}

library(dplyr)
library(foreach)
dtb <- NULL
dtb_base <- NULL


for (batch in 0:batch_length) {
  print(sprintf("%s, %s", Sys.time(), batch))
  store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
  data <- readRDS(store_path)
  data_tibble <- nl2tibble(data$nl)
  data_tibble <- hts.eval2(data_tibble, metrics, data$tts, data$bts)

  data_tibble <- data_tibble %>% select(-rf) %>% mutate(batch = batch)
  dtb <- rbind(dtb, data_tibble)
  dtb_base <- rbind(dtb_base, as_tibble(hts.evalbase(data, metrics)) %>%
                      mutate(batch = batch))
}

saveRDS(list(base = dtb_base, dtb = dtb), sprintf("%s/%s/eval.rds", path, bfmethod))



