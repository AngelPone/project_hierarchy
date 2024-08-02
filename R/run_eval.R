args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- "ets"

dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
time_length <- NROW(dt$data)
forecast_horizon <- 12
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

source("R/metrics.R")
source("R/baseforecast.R")


hts.eval2 <- function(df, tts, bts, S, combination, base, comb_permute=NULL) {
  tts <- tts %*% t(S)
  bts <- bts %*% t(S)
  
  df <- add_row(df, representor="", cluster="combination1", 
                distance="", S=NULL, 
                rf=list(list(mint=combination[[1]])), other=NULL) 
  
  df <- add_row(df, representor = "", cluster="combination2", S=NULL, 
                distance = "",
                rf=list(list(mint=combination[[2]])), other=NULL)
  
  if (!is.null(comb_permute)) {
    for (i in 1:100) {
      df <- add_row(df, representor="", cluster=paste0("permute-combination1-", i),
                    distance="", S=NULL, rf=list(list(mint=comb_permute[[i]])), other=NULL)
    }
  }
  
  
  df[["rmsse"]] <- lapply(iterators::iter(df$rf), function(g) {
      c <- g[['mint']][, 2:(m+1),drop=FALSE]
      c <- c %*% t(S)
      sapply(1:NCOL(c), function(x) { metric.rmsse(tts[,x], c[,x], bts[,x]) })
  })
  
  
  basef_middle <- furrr::future_map(as.list(iterators::iter(bts[,2:(n-m)], by="column")),
                                    \(x) as.numeric(f.ets(x, 12, 12)$basef))
  basef_middle <- do.call(cbind, basef_middle)
  basef <- cbind(base[,1], basef_middle, base[,2:(m+1)])
  basef_rmsse <- sapply(1:NCOL(bts),
                             function(x) {
                               metric.rmsse(tts[,x], basef[,x], bts[,x])
                             } )
  
  df <- add_row(df, representor="", cluster="base", distance="",
                S=NULL, rf=NULL, rmsse=list(basef_rmsse))
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


library(dplyr)
library(foreach)
library(furrr)
plan(multisession(workers = 8))
dtb <- NULL




combination <- readRDS(sprintf("%s/ets/combination.rds", path))
if (path == "mortality") {
  combination_permute <- readRDS("mortality/ets/combination_permute.rds")
}
for (batch in 0:batch_length) {
  print(sprintf("%s, %s", Sys.time(), batch))
  store_path <- sprintf("%s/ets/batch_%s.rds", path, batch)
  data <- readRDS(store_path)
  data_tibble <- nl2tibble(data$nl)
  c_p <- NULL
  if (path == "mortality") {
    c_p <- combination_permute[[batch+1]]
  }
  data_tibble <- hts.eval2(data_tibble, data$tts, data$bts, dt$S, combination[[batch+1]], data$basef, c_p)

  data_tibble <- data_tibble %>% select(-rf, -other) %>% mutate(batch = batch)
  dtb <- rbind(dtb, data_tibble)
}

saveRDS(dtb, sprintf("%s/ets/eval.rds", path))



