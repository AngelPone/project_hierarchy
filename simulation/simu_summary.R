rm(list = ls())
dt <- readRDS("simulation/simulation.rds")
library(dplyr)
source("simulation/utils.R")

forecast_horizon <- 1


metric.rmsse <- function(obs, pred, hist) {
  sqrt(mean((obs - pred)^2) / mean(diff(hist, 12)^2))
}

rmsse <- function(pred, obs, hist) {
  obs <- cbind(rowSums(obs), obs)
  hist <- cbind(rowSums(hist), hist)
  mean(sapply(1:NCOL(pred), function(x) {
    metric.rmsse(
      pred[1:forecast_horizon, x, drop = FALSE],
      obs[1:forecast_horizon, x, drop = FALSE],
      hist[, x, drop = FALSE]
    )
  }))
}


no_hierarchy <- 1
best_hierarchy <- 2
permute_best <- 3:102

permute_trend <- 103:202
trend_dir <- 303

season <- 304
permute_season <- 203:302

trend_exis <- 305
permute_trend_exis <- 306:405

evaluate_idx <- function(idx) {
  if (length(idx) > 1) {
    output <- do.call(cbind, lapply(idx, function(x) {
      evaluate_idx(x)
    }))
    return(output)
  }
  sapply(1:500, function(x) {
    rmsse(
      dt$acc[[x]][[idx]],
      t(dt$series[[x]][, 133:144]),
      t(dt$series[[x]][, 1:132])
    )
  })
}

compare_random <- function(idx_orig, idx_random) {
  do.call(cbind, lapply(c(idx_orig, idx_random), evaluate))
}

test <- function(mat, name) {
  colnames(mat) <- c(name, 1:100)
  nemenyi(mat, plottype = "vmcb", target = name)
}

# natural hierarchy vs its counterpart
jpeg("figures/Figure11.jpg", width = 600 *1.4, height = 300*1.4)
par(mar=c(4,18,3,2))
natural_ <- evaluate_idx(best_hierarchy)
natural_test <- test(cbind(natural_, evaluate_idx(permute_best)), "Cluster-trend-season")
dev.off()


season_ <- evaluate_idx(season)
trend_dir_ <- evaluate_idx(trend_dir)
trend_exis_ <- evaluate_idx(trend_exis)

two_level <- evaluate_idx(1)


calculate_base <- function() {
  if (file.exists("simulation/base.rds")) {
    return(sapply(readRDS("simulation/base.rds"), function(x) {
      x
    }))
  }
  library(forecast)
  base_ <- lapply(1:500, function(x) {
    all_series <- t(rbind(colSums(dt$series[[x]]), dt$series[[x]]))
    train <- all_series[1:132, ]
    test <- all_series[133:(133 + forecast_horizon), , drop = FALSE]
    fcasts <- lapply(iterators::iter(train, by = "column"), function(y) {
      mdl <- ets(ts(y, frequency = 12))
      as.numeric(
        forecast(mdl, h = forecast_horizon)$mean
      )
    }) %>% do.call(cbind, .)
    rmsse(fcasts, test, train)
  })
  saveRDS(base_, "simulation/base.rds")
  return(sapply(base_, function(x) {
    x
  }))
}


print("evaluating base ...")
base_ <- calculate_base()


cluster_mat <- cbind(base_, two_level, natural_, trend_dir_, trend_exis_, season_)
colnames(cluster_mat) <- c("Base", "Two-level", "Cluster-trend-season", "Cluster-trend1", "Cluster-trend2", "Cluster-season")
clusters_tbl <- data.frame(
  method = c("Base", "Two-level", "Cluster-trend-season", "Cluster-trend1", "Cluster-trend2", "Cluster-season"),
  rmsse = colMeans(cluster_mat) * 100
) %>%
  arrange(rmsse) %>%
  mutate(rmsse = round(rmsse, digits = 2)) %>%
  write.csv("figures/Table7.csv")

  

jpeg("figures/Figure10.jpg", width = 600, height = 300)
tsutils::nemenyi(cluster_mat, plottype = "vmcb")
dev.off()
 




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


library(ggplot2)

visualizeSimPCA(dt$series[[1]])
ggsave("figures/Figure8.jpg", width = 6, height = 3)



jpeg("figures/Figure9.jpg", width = 1200, height = 600)
visualizeSimSeries(dt$series[[1]])
dev.off()
 
