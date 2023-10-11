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
  ts <- get(paste0("representator.", representor))(dataset$data)
  pca_res <- princomp(cor(ts))$loadings[,1:2]
  
  cluster <- dataset$output$cluster == cluster
  distance <- dataset$output$distance == distance
  representor <- dataset$output$representator == representor
  
  cluster_S <- dataset$output$other[[which(cluster & distance & representor)]]$S
  cluster_S <- as.matrix(cluster_S)
  
  grp <- vector("numeric", NROW(pca_res))
  
  for (i in 1:NROW(cluster_S)) {
    grp[which(cluster_S[i,] == 1)] <- i
  }
  pca_res <- as_tibble(pca_res) %>% mutate(group = as.character(grp))
  
  
  ggplot(pca_res) +
    geom_point(aes(x = Comp.1, y=Comp.2, color=group, group=group))
}


dmrvisualize(tourism, "ts", "kmedoids-5", "euclidean")
dmrvisualize(mortality, "ts", "kmedoids-5", "euclidean")


dmrvisualize(mortality, "error", "kmedoids-5", "negcor")
dmrvisualize(tourism, "error", "kmedoids-5", "negcor")

dmrvisualize(mortality, "ts", "kmedoids-5", "euclidean")



