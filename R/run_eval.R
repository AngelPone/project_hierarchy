args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

metrics <- c("rmse", "mae", "rmsse", "mase")
source("R/metrics.R")

hts.eval <- function(df, metrics, tts, bts) {
  tts <- cbind(rowSums(tts), tts)
  bts <- cbind(rowSums(bts), bts)
  
  # random
  randoms <- unique(df$cluster[which(startsWith(df$cluster, "random"))])
  randoms <- c(randoms, "hcluster-random")
  df_random <- list()
  df_random$S <- vector("list", 3*length(randoms) + 1)
  df_random$other <- vector("list", 3*length(randoms) + 1)
  df_random$cluster <- c()
  df_random$rf <- list()
  df_random$representor <- rep("", 3*length(randoms) + 1)
  df_random$distance <- rep("", 3*length(randoms) + 1)
  for (rd in randoms) {
    tmpdt <- df %>% filter(cluster == rd)
    for (random_n in c(10, 20, 50)) {
      avg_rf <- list()
      for (rf_method in c("ols", "wlss", "wlsv", "mint")) {
        avg_rf_method <- 
          do.call(abind::abind, list(lapply(tmpdt$rf, function(x) x[[rf_method]]), along=0))
        avg_rf[[rf_method]] <- 
          apply(avg_rf_method[1:random_n, ,], c(2, 3), mean)
      }
      df_random$cluster <- c(df_random$cluster, paste0(rd, "-", random_n))
      df_random$rf <- append(df_random$rf, list(avg_rf))
    }
  }
  # cluster average
  df_random$cluster <- c(df_random$cluster, "cluster-average")
  avg_rf <- list()
  tmpdt <- df %>% filter(cluster %in% c("Kmedoids-dr", "hcluster-dr")) %>%
    rbind(df %>% filter(distance == "dtw", representor %in% c("ts", "error")))
  for (rf_method in c("ols", "wlss", "wlsv", "mint")) {
    avg_rf_method <- 
      do.call(abind::abind, list(lapply(tmpdt$rf, function(x) x[[rf_method]]), along=0))
    avg_rf[[rf_method]] <- 
      apply(avg_rf_method, c(2, 3), mean)
  }
  df_random$rf <- append(df_random$rf, list(avg_rf))
  
  df <- df %>% filter(!startsWith(cluster, "random")) %>%
    filter(cluster != "hcluster-random") %>%
    rbind(as_tibble(df_random))
  
  for (metric in metrics) {
    accuracy_method <- get(paste0("metric.", metric))
    
    df[[metric]] <- foreach::foreach(g=iterators::iter(df$rf)) %dopar%  {
      lapply(g, function(c) { 
        sapply(1:NCOL(c), function(x) { accuracy_method(tts[,x], c[,x], bts[,x]) } )
      })
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
  tts <- cbind(rowSums(dt$tts), dt$tts)
  bts <- cbind(rowSums(dt$bts), dt$bts)
  output <- list()
  for (metric in metrics) {
    accuracy_method <- get(paste0("metric.", metric))
    output[[metric]] <-
      list(sapply(1:NCOL(bts), function(x) { accuracy_method(tts[,x], dt$basef[,x], bts[,x]) } ))
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
  data_tibble <- hts.eval(data_tibble, metrics, data$tts, data$bts)

  data_tibble <- data_tibble %>% select(-rf) %>% mutate(batch = batch)
  dtb <- rbind(dtb, data_tibble)
  dtb_base <- rbind(dtb_base, as_tibble(hts.evalbase(data, metrics)) %>%
                      mutate(batch = batch))
}

saveRDS(list(base = dtb_base, dtb = dtb), sprintf("%s/%s/eval.rds", path, bfmethod))



