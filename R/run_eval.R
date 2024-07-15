args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- "ets"

# cl <- parallel::makeCluster(8)
# doParallel::registerDoParallel(cl)

dt <- readRDS(sprintf("%s/data.rds", path))
n <- NROW(dt$S)
m <- NCOL(dt$S)
print(sprintf("%s dataset has %s series and %s bottom series", path, n, m))
time_length <- NROW(dt$data)
forecast_horizon <- as.integer(args[[2]])
frequency <- 12
batch_length <- time_length - 96 - forecast_horizon

metrics <- c("rmsse")
source("R/metrics.R")


hts.eval2 <- function(df, metrics, tts, bts, S) {
  if (is.null(dim(tts))) {
    tts <- S %*% tts
    tts <- matrix(tts, nrow = 1)
  } else {
    tts <- (tts %*% t(S))[1:forecast_horizon,,drop=FALSE]
  }

  bts <- bts %*% t(S)

  data_tibble <- df %>% filter(representor != "")

  data_tibble$permute <- sapply(data_tibble$cluster, function(x){
    if (!startsWith(x, "permute")) {
      return (0)
    }
    x <- strsplit(x, "-")[[1]]
    x <- x[length(x)]
    as.integer(x)
  })
  data_tibble$cluster <- sapply(data_tibble$cluster, function(x){
    if (!startsWith(x, "permute")) {
      return (x)
    }
    split_x <- strsplit(x, "-")[[1]]
    x <- stringi::stri_replace_all_fixed(x, "permute-", "")
    x <- stringi::stri_replace_all_fixed(x, paste0("-", split_x[length(split_x)]), "")
    x
  })

  # compute combination of reconciled forecasts
  avg <- data_tibble %>% select(representor, cluster, distance, rf, permute) %>%
    arrange(permute) %>%
    group_by(permute) %>%
    tidyr::nest(rf = -"permute") %>%
    mutate_at("rf", purrr::map, function(x) {
      if (dim(x)[1] != 12) {
        return(NULL)
      }
      output <- lapply(c("ols", "wlss", "wlsv", "mint"), function(m) {
        sapply(x$rf, function(g) { g[[m]][1:forecast_horizon,,drop=FALSE] }, simplify = "array") %>%
          apply(c(1,2), mean)
      })
      names(output) <- c("ols", "wlss", "wlsv", "mint")
      output
    }) %>% rowwise() %>%
    filter(!is.null(rf)) %>% ungroup() %>%
    mutate(permute = paste0(ifelse(permute > 0, "permute-", ""),
                            "average",
                            ifelse(permute > 0, paste0("-", permute), ""))) %>%
    mutate(representor="", distance="", other=list(NULL), S=list(NULL)) %>%
    rename(cluster=permute)

  df <- df %>% rbind(avg)


  for (metric in metrics) {
    accuracy_method <- get(paste0("metric.", metric))

    df[[metric]] <- foreach::foreach(g=iterators::iter(df$rf)) %do%  {
      c <- g[['mint']][1:forecast_horizon, 2:(m+1),drop=FALSE]
      c <- c %*% t(S)
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

hts.evalbase <- function(dt, metrics, S) {
  if (is.null(dim(dt$tts))) {
    tts <- S %*% tts
    tts <- matrix(tts, nrow = 1)
  } else {
    tts <- (tts %*% t(S))[1:forecast_horizon,,drop=FALSE]
  }

  bts <- bts %*% t(S)
  output <- list()
  basef_middle <- foreach(iterators::iter(x=bts[,2:(n-m)], by="column")) %dopar% {
    x <- ts(x, frequency = 12)
    mdl <- forecast::ets(x)
    fcasts <- forecast::forecast(mdl, h=12)
    as.numeric(fcasts$mean)
  } %>% do.call(cbind)
  basef <- cbind(dt$basef[,1], basef_middle, dt$basef[,2:(m+1)])
  for (metric in metrics) {
    accuracy_method <- get(paste0("metric.", metric))
    output[[metric]] <-
      list(sapply(1:NCOL(bts),
                  function(x) {
                    accuracy_method(tts[,x], basef[1:forecast_horizon,x], bts[,x])
                    } ))
  }
  output
}

library(dplyr)
library(foreach)
dtb <- NULL
dtb_base <- NULL


for (batch in 0:batch_length) {
  print(sprintf("%s, %s", Sys.time(), batch))
  store_path <- sprintf("%s/ets/batch_%s.rds", path, batch)
  data <- readRDS(store_path)
  data_tibble <- nl2tibble(data$nl)
  data_tibble <- hts.eval2(data_tibble, metrics, data$tts, data$bts, dt$S)

  data_tibble <- data_tibble %>% select(-rf) %>% mutate(batch = batch)
  dtb <- rbind(dtb, data_tibble)
  dtb_base <- rbind(dtb_base, as_tibble(hts.evalbase(data, metrics)) %>%
                      mutate(batch = batch))
}

saveRDS(list(base = dtb_base, dtb = dtb), sprintf("%s/%s/eval_%s.rds", path, bfmethod, forecast_horizon))



