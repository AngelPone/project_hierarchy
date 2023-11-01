library(dplyr)
library(tidyr)
dt <- readRDS("simulation/simulation.rds")

acc <- dt$acc %>%
  nest(values = c("batch", "rmse")) %>%
  mutate_at("values", purrr::map, function(x) {
    x <- arrange(x, batch)
    do.call(c, x$rmse)
  })

acc_mat <- do.call(cbind, acc$values) %>%
  `colnames<-`(acc$methods)



FREQUENCIES <- c(12, 4, 12, 4)
visualizeSimSeries <- function(series) {
  par(mfrow = c(2,2))
  grpidx <- c("1", "2", "3", "4")
  for (i in 0:3) {
    plot(ts(series[20*i+1,], frequency = FREQUENCIES[i+1]), ylab = "series", main = paste0("group", grpidx[i+1]))
  }
  par(mfrow = c(1,1))
}

visualizeSimPCA <- function(series) {
  pca <- prcomp(series, scale.=TRUE)
  grpidx <- c("1", "2", "3", "4")
  grps <- rep(grpidx, each=20)
  ggplot(mapping = aes(x = pca$x[,1], y=pca$x[,2], color = grps, group=grps)) +
    geom_point() +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") +
    labs(color = "groups")
}

pdf("simulation/mcb.pdf", width = 10, height = 8)
tsutils::nemenyi(acc_mat, plottype = "vmcb")
dev.off()

pdf("simulation/PCA.pdf", width = 8, height = 8)
visualizeSimPCA(dt$series[[1]])
dev.off()


pdf("simulation/example.pdf", width = 16, height = 8)
visualizeSimSeries(dt$series[[1]])
dev.off()
