rm(list = ls())
dt <- readRDS("simulation/simulation.rds")
library(dplyr)
source("R/utils.R")

forecast_horizon <- 1

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
pdf("manuscript/figures/simulation/P3_c_vs_pc.pdf", width = 10, height = 5)
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

calculate_comb <- function() {
  if (file.exists("simulation/comb.rds")) {
    return(readRDS("simulation/comb.rds"))
  }
  method_idx <- c(2, 303, 304, 305)
  all_rmsse <- NULL
  for (i in seq_along(dt$acc)) {
    print(i)
    tts <- t(dt$series[[i]][, 133:144])
    bts <- t(dt$series[[i]][, 1:132])
    f_orig_ <- apply(simplify2array(dt$acc[[i]][method_idx]), c(1, 2), mean)
    f_permu_ <- lapply(1:100, function(j) {
      apply(simplify2array(
        dt$acc[[i]][c(permute_best[j], permute_trend[j], permute_season[j], permute_trend_exis[j])]
      ), c(1, 2), mean)
    })

    rmsse_orig <- rmsse(f_orig_, tts, bts)
    rmsse_permu_ <- sapply(iterators::iter(f_permu_), function(g) {
      rmsse(g, tts, bts)
    })
    all_rmsse <- rbind(all_rmsse, c(rmsse_orig, rmsse_permu_))
  }
  saveRDS(all_rmsse, "simulation/comb.rds")
  all_rmsse
}

print("Evaluating comb ...")
comb_rmsse <- calculate_comb()
# # 51
pdf("manuscript/figures/simulation/P4.pdf", width = 8, height = 6)
comb_test <- test(comb_rmsse, "Combination")
dev.off()


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
  write.csv("manuscript/figures/simulation/simulation_methods.csv")

  

pdf("manuscript/figures/simulation/P3_mcb.pdf", width = 8, height = 4)
tsutils::nemenyi(cluster_mat, plottype = "vmcb")
dev.off()

output <- list()
output$cluster <- c(
  which(names(natural_test$means) == "Cluster-trend-season"),
  mean(natural_)
)

output$comb <- c(
  which(names(comb_test$means) == "Combination"),
  mean(comb_rmsse[, 1])
)
data.frame(output) %>% write.csv("manuscript/figures/simulation_ranktbl.csv")
