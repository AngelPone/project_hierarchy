batch <- 11
library(dplyr)
library(ggplot2)
mortality <- readRDS(sprintf("mortality/store_%s.rds", batch))$data
tourism <- readRDS(sprintf("tourism/arima/store_%s.rds", batch))$data

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
visual(tourism, sample(304, 4))
visual(mortality, sample(98, 4))



# dimensional reduction visualization
tourism <- readRDS("tourism/arima/store_0.rds")
mortality <- readRDS(sprintf("mortality/arima/store_%s.rds", batch))
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
  plot(ts(dataset$data$bts[, idx+1], frequency = 12), main = sprintf("time series entropy = %.2f", entropy(dataset$data$bts[, idx+1], order = order)))
  acf(dataset$data$resid[, idx], main = sprintf("acf of residual entropy = %.2f", entropy(dataset$data$resid[,idx], order = order)))
  acf(dataset$data$bts[, idx], main = "acf of time series")
}

visualizeRSA(tourism, sample(98, 1))


visualizeRSA(tourism, 33, order = 12)
visualizeRSA(tourism, 78, order = 12)
visualizeRSA(tourism, 93, order = NULL)
visualizeRSA(tourism, 79, order = 12)




visualizeGroupCOR <- function(dataset, representor, distance, cluster) {
  cluster_ <- dataset$output$cluster == cluster
  representor_ <- dataset$output$representator == representor
  distance_ <- dataset$output$distance == distance
  cluster_S <- as.matrix(dataset$output$other[[which(cluster_ & representor_ & distance_)]]$S)
  
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


visualizeGroupCOR(mortality, "error", "uncorrelation", "kmedoids-4")
dmrvisualize(mortality, "error", "kmedoids-4", "euclidean")




