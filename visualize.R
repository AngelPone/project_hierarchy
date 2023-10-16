batch <- 11
library(dplyr)
library(ggplot2)
bfmethod <- "ets"

mortality <- readRDS(sprintf("mortality/%s/store_%s.rds", bfmethod, batch))
tourism <- readRDS(sprintf("tourism/%s/store_%s.rds", bfmethod, batch))

mortality_arima <- readRDS(sprintf("mortality/%s/store_%s.rds", "arima", batch))
tourism_arima <- readRDS(sprintf("tourism/%s/store_%s.rds", "arima", batch))


visual <- function(dt, idx) {
  
  resid <- unname(as.matrix(dt$resid[,idx+1]))
  ts <- unname(as.matrix(dt$bts[,idx]))
  
  resid <- apply(resid, 2, function(x){
    (x - mean(x))/sd(x)
  })
  ts <- apply(ts, 2, function(x){
    (x-mean(x))/sd(x)
  })
  
  colnames(ts) <- paste("Series", 1:NCOL(ts))
  colnames(resid) <- paste("Series", 1:NCOL(ts))
  
  data.frame(resid) %>% mutate(type = "error", x=1:NROW(resid)) %>%
    rbind(data.frame(ts) %>% mutate(type = "ts", x=1:NROW(resid))) %>%
    tidyr::pivot_longer(cols = -c("type", "x")) %>%
    ggplot(mapping = aes(x=x, y=value)) +
    geom_line() +
    facet_wrap(~ type + name, nrow = 2)
}
visual(tourism$data, sample(304, 4))
visual(mortality$data, sample(98, 4))



# dimensional reduction visualization
dmrvisualize <- function(dataset, representor, cluster, distance) {
  dataset$data$features <- dataset$features
  ts <- get(paste0("representator.", representor))(dataset$data)
  pca_res <- princomp(cor(ts))$loadings[,1:2]
  
  cluster_ <- dataset$output$cluster == cluster
  distance_ <- dataset$output$distance == distance
  representor_ <- dataset$output$representator == representor
  
  cluster_S <- dataset$output$other[[which(cluster_ & distance_ & representor_)]]$S
  cluster_S <- as.matrix(cluster_S)
  
  grp <- vector("numeric", NROW(pca_res))
  
  for (i in 1:NROW(cluster_S)) {
    grp[which(cluster_S[i,] == 1)] <- i
  }
  pca_res <- as_tibble(pca_res) %>% mutate(group = as.character(grp))
  
  
  ggplot(pca_res) +
    geom_point(aes(x = Comp.1, y=Comp.2, color=group, group=group)) +
    ggtitle(paste0(representor, " ", cluster, " ", distance))
}


gridExtra::grid.arrange(
  dmrvisualize(tourism, "ts", "kmedoids-5", "euclidean"),
  dmrvisualize(tourism, "ts", "kmedoids-5", "dtw"),
  dmrvisualize(tourism, "error", "kmedoids-5", "euclidean"),
  dmrvisualize(tourism, "error", "kmedoids-5", "dtw"),
  dmrvisualize(tourism, "ts.features", "kmedoids-5", "euclidean"),
  dmrvisualize(tourism, "ts.features", "kmedoids-5", "dtw"),
  dmrvisualize(tourism, "error.features", "kmedoids-5", "euclidean"),
  dmrvisualize(tourism, "error.features", "kmedoids-5", "dtw")
)

gridExtra::grid.arrange(
  dmrvisualize(mortality, "ts", "kmedoids-5", "euclidean"),
  dmrvisualize(mortality, "ts", "kmedoids-5", "dtw"),
  dmrvisualize(mortality, "error", "kmedoids-5", "euclidean"),
  dmrvisualize(mortality, "error", "kmedoids-5", "dtw"),
  dmrvisualize(mortality, "ts.features", "kmedoids-5", "euclidean"),
  dmrvisualize(mortality, "ts.features", "kmedoids-5", "dtw"),
  dmrvisualize(mortality, "error.features", "kmedoids-5", "euclidean"),
  dmrvisualize(mortality, "error.features", "kmedoids-5", "dtw")
)



# visualize Residual, series and acf
visualizeRSA <- function(dataset, idx, order = NULL) {
  par(mfrow = c(2,2))
  plot(dataset$data$resid[, idx+1], main = sprintf("residual index=%s", idx))
  plot(ts(dataset$data$bts[, idx], frequency = 12), main = sprintf("time series entropy = %.2f,unitroot = %.2f", entropy(dataset$data$bts[, idx], order = order), tsfeatures::unitroot_kpss(dataset$data$bts[, idx])))
  acf(dataset$data$resid[, idx+1], main = sprintf("acf of residual entropy = %.2f", entropy(dataset$data$resid[,idx], order = order)))
  acf(dataset$data$bts[, idx], main = "acf of time series")
}


visualizeRSA(tourism, sample(304, 1))


visualizeRSA(tourism, 72, order = 12)
visualizeRSA(tourism, 80, order = 12)
visualizeRSA(tourism, 1, order = NULL)
visualizeRSA(tourism, 5, order = NULL)
visualizeRSA(mortality, sample(98, 1), order = NULL)


visualizeGroupCOR <- function(dataset, representor, distance, cluster, idx = NULL) {
  
  if (!is.null(idx)) {
    cluster_S <- as.matrix(dataset$output$other[[idx]]$S)
  } else {
    cluster_ <- dataset$output$cluster == cluster
    representor_ <- dataset$output$representator == representor
    distance_ <- dataset$output$distance == distance
    cluster_S <- as.matrix(dataset$output$other[[which(cluster_ & representor_ & distance_)]]$S)
  }
  
  par(mfrow = c(1, 2))
  for (i in 1:NROW(cluster_S)) {
    total <- dataset$data$bts %*% matrix(cluster_S[i,], ncol = 1)
    series <- cbind(total, dataset$data$bts[,which(cluster_S[i, ] == 1)])
    resid <- dataset$data$resid[,which(cluster_S[i, ] == 1) + 1]
    colnames(resid) <- which(cluster_S[i, ] == 1)
    colnames(series) <- c("total", which(cluster_S[i, ] == 1))
    corrplot::corrplot(cor(series), method = "color")
    corrplot::corrplot(cor(resid), method = "color")
  }
}

accuracy_order <- function(dataset, rf_method, accuracy_method){
  lapply(dataset$output$accuracy[which(!startsWith(mortality$output$cluster, "random-single"))[-1]], function(x){
    c(x[[accuracy_method]][[rf_method]]$total, 
      x[[accuracy_method]][[rf_method]]$bottom)
  }) %>% do.call(rbind, .) %>%
    rbind(c(dataset$output$accuracy[[1]][[accuracy_method]]$base$total,
            dataset$output$accuracy[[1]][[accuracy_method]]$base$bottom),
          .)
}

method_comb <- function(dataset, idx) {
  idx <- which(!startsWith(dataset$output$cluster, "random-single"))[idx]
  list(representor = dataset$output$representator[idx],
    distance = dataset$output$distance[idx],
    cluster = dataset$output$cluster[idx])
}




visualizeGroupCOR(mortality, best_$representor, best_$distance, best_$cluster, idx = NULL)

visualizeSeriesError <- function(dataset, order = NULL, method = NULL) {
  if (is.null(method)) {
    acc_rank <- accuracy_order(dataset, "mint", "rmse")
    rank_rmse_mint <- cbind(acc_rank[,1], 
                            acc_rank[,2:NCOL(acc_rank)]) %>% t() %>%
      apply(1, rank) %>%
      rowMeans()
    method <- method_comb(dataset, which(order(rank_rmse_mint) == order))
  }
  print(sapply(method, function(x){x}))
  cluster_ <- dataset$output$cluster == method$cluster
  representor_ <- dataset$output$representator == method$representor
  distance_ <- dataset$output$distance == method$distance
  cluster_S <- as.matrix(dataset$output$other[[which(cluster_ & representor_ & distance_)]]$S)
  par(mfrow = c(1,2))
  for (i in 1:NROW(cluster_S)) {
    if (sum(cluster_S[i,]) <= 10) {
      series <- ts(dataset$data$bts[,which(cluster_S[i,] == 1)], frequency = 12)
      # colnames(series) <- colnames(mortality_S)[which(cluster_S[i,] == 1)]
      resids <- dataset$data$resid[,which(cluster_S[i,] == 1) + 1]
      # colnames(resids) <- colnames(mortality_S)[which(cluster_S[i,] == 1)]
      plot(series, main = sprintf("series Group = %s", i))
      plot(resids, main = sprintf("residuals Group = %s", i))
    }
  }
}

# mortality_S <- read.csv("data/S.csv")
# colnames(mortality_S) <- stringi::stri_replace_all_fixed(colnames(mortality_S), ".", "-")
# mortality_S <- mortality_S[,2:NCOL(mortality_S)]


visualizeSeriesError(mortality, 2)
visualizeSeriesError(mortality_arima, method = list(representor = "ts.features",
                                        distance = "dtw",
                                        cluster = "kmedoids-10"))

visualizeSeriesError(tourism, 3)
visualizeSeriesError(tourism_arima, 6)
