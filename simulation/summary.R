library(dplyr)
library(tidyr)
library(ggplot2)

setting <- 2

dt <- readRDS(sprintf("simulation/simulation%s.rds", setting))

acc <- dt$acc %>%
  nest(values = c("batch", "rmse")) %>%
  mutate_at("values", purrr::map, function(x) {
    x <- arrange(x, batch)
    do.call(c, x$rmse)
  })

acc$methods <- c("C0", "C1", "A2", "random-20", "random-50", "C2", "C4", "C3", "A1", "Base")
acc <- acc %>% filter(!startsWith(methods, "random"))

acc_mat <- do.call(cbind, acc$values) %>%
  `colnames<-`(acc$methods)




visualizeSimSeries <- function(series) {
  par(mfrow = c(3,2))
  grpidx <- c("1", "2", "3", "4", "5", "6")
  for (i in 0:5) {
    plot(ts(series[20*i+1,], frequency = 12), ylab = "series", main = paste0("Cluster", grpidx[i+1]))
  }
  par(mfrow = c(1,1))
}

visualizeSimPCA <- function(series) {
  pca <- prcomp(series, scale.=TRUE)
  grpidx <- c("1", "2", "3", "4", "5", "6")
  grps <- rep(grpidx, each=20)
  ggplot(mapping = aes(x = pca$x[,1], y=pca$x[,2], color = grps, group=grps)) +
    geom_point() +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") +
    labs(color = "Clusters")
}

pdf(sprintf("manuscript/figures/simu_mcb%s.pdf", setting), width = 5, height = 5)
tsutils::nemenyi(acc_mat, plottype = "vmcb")
dev.off()

pdf("manuscript/figures/simu_pca.pdf", width = 4, height = 2.5)
visualizeSimPCA(dt$series[[1]])
dev.off()


pdf("manuscript/figures/simu_example.pdf", width = 16, height = 8)
visualizeSimSeries(dt$series[[1]])
dev.off()
