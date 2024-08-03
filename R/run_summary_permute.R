source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)
library(tsutils)

source("R/metrics.R")
source("R/expr_utils.R")




# natural vs its randomization
P2 <- function(dt, path) {

  # MCB Test based on hierarchy rmsse
  pdf(sprintf("manuscript/figures/%s/natural_vs_pn.pdf", path), width = 8, height = 6)
  par(mar=c(4,14,3,2))
  natural_hierarchy <- mcb_hierarchy_rmsse(
    dt %>% filter(cluster == "natural"),
    dt %>% filter(startsWith(cluster, "permute-natural")),
    "Natural"
  )
  dev.off()
  rank <- which(names(natural_hierarchy$means) == "Natural")
  itl1 <- natural_hierarchy$means - natural_hierarchy$cd / 2
  itl2 <- natural_hierarchy$means + natural_hierarchy$cd / 2
  sig_better <-
    length(which(itl1 > itl2["Natural"]))
  sig_worse <-
    length(which(itl2 < itl1["Natural"]))

  write(sprintf("Natural ranks %s in its 100 twins, significantly better than %s, significantly worse than %s",
                rank, sig_better, sig_worse), 
        sprintf("manuscript/figures/%s/natural_vs_pn.txt", path))
}


P3 <- function(dt, path) {
  bench_rmsse <- dt %>%
    filter(!startsWith(cluster, "permute"), !startsWith(cluster, "combination")) %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse)

  pdf(sprintf("manuscript/figures/%s/mcb_benchmarks.pdf", path), 8, 6)
  par(mex = 1.1)
  bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()

  best_ <- bench_rmsse %>%
    group_by(representor, distance, cluster) %>%
    summarise(rmsse = mean(rmsse), .groups = "drop") %>%
    filter(representor != "") %>%
    arrange(rmsse)

  best_ <- best_[1, ]

  best_name <- method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])


  pdf(sprintf("manuscript/figures/%s/cluster_vs_pc.pdf", path), 8, 6)
  par(mar=c(4,14,3,2))
  cluster_hierarchy_rmsse <- mcb_hierarchy_rmsse(
    dt %>% filter(
      representor == best_$representor[[1]],
      cluster == best_$cluster[[1]],
      distance == best_$distance[[1]]
    ),
    dt %>% filter(
      representor == best_$representor[[1]],
      startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
      distance == best_$distance[[1]]
    ),
    method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])
  )
  dev.off()


  # calculate rmsse
  bench_rmsse <- dt %>%
    filter(!startsWith(cluster, "permute"), !startsWith(cluster, "combination")) %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    mutate(top=rmsse[1], middle=mean(rmsse[2:(n-m)]), bottom=mean(rmsse[(n-m+1):n])) %>% 
    ungroup() %>%
    group_by(method) %>%
    summarise(across(c(top, middle, bottom), \(x) round(mean(x), digits = 4)))

  rank <- which(names(cluster_hierarchy_rmsse$means) == best_name)
  itl1 <- cluster_hierarchy_rmsse$means - cluster_hierarchy_rmsse$cd / 2
  itl2 <- cluster_hierarchy_rmsse$means + cluster_hierarchy_rmsse$cd / 2
  sig_better <- length(which(itl1 > itl2[best_name]))
  sig_worse <- length(which(itl2 < itl1[best_name]))
  
  write(sprintf("The best cluster ranks %s in its 100 twins, significantly better than %s, significantly worse than %s",
                rank, sig_better, sig_worse), 
        sprintf("manuscript/figures/%s/cluster_vs_pc.txt", path))
  
  methods <- c(
    "Base", "Two-level", "Natural", 
    "TS-EUC-ME", "ER-EUC-ME", "TSF-EUC-ME", "ERF-EUC-ME",
    "TS-EUC-HC", "ER-EUC-HC", "TSF-EUC-HC", "ERF-EUC-HC",
    "TS-DTW-ME", "TS-DTW-HC", "ER-DTW-ME", "ER-DTW-HC")
    
  bench_rmsse[match(methods, bench_rmsse$method),] %>%
    write.csv(sprintf("manuscript/figures/%s/rmsse_benchmarks.csv", path))
}

P4 <- function(dt, path) {
  bench_rmsse <- dt %>%
    filter(!startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse), method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse)

  pdf(sprintf("manuscript/figures/%s/mcb_combination.pdf", path), 6, 8)
  bench_rmsse %>%
    tidyr::pivot_wider(names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()


  methods <- c(
    "Base", "Two-level", "Natural", 
    "TS-EUC-ME", "ER-EUC-ME", "TSF-EUC-ME", "ERF-EUC-ME",
    "TS-EUC-HC", "ER-EUC-HC", "TSF-EUC-HC", "ERF-EUC-HC",
    "TS-DTW-ME", "TS-DTW-HC", "ER-DTW-ME", "ER-DTW-HC",
     "Combination", "Stack"
  )
  bench_rmsse <- dt %>%
    filter(!startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(top = rmsse[1], middle=mean(rmsse[2:(n-m)]), bottom=mean(rmsse[(n-m+1):n]),
           method = method_name(representor, distance, cluster)) %>%
    ungroup() %>%
    select(method, batch, top, middle, bottom) %>%
    group_by(method) %>%
    summarise(across(c(top, middle, bottom), \(x) round(mean(x), 4)))
  write.csv(
    bench_rmsse[match(methods, bench_rmsse$method),],
    sprintf("manuscript/figures/%s/combination.csv", path))
}

P5 <- function(dt, path) {
  
  pdf(sprintf("manuscript/figures/%s/comb_vs_pc.pdf", path), 6, 8)
  mcb <- mcb_hierarchy_rmsse(
    dt %>% filter(cluster == "combination1"),
    dt %>% filter(startsWith(cluster, "permute-combination1")),
    "Combination"
  )
  dev.off()
}


for (path in c("tourism", "mortality")) {
  dir.create(paste0("manuscript/figures/", path))
  dt <- readRDS(sprintf("%s/ets/eval.rds", path))
  orig_data <- readRDS(sprintf("%s/data.rds", path))
  n <- NROW(orig_data$S)
  m <- NCOL(orig_data$S)
  print("Part 2...")
  P2(dt, path)
  print("Part 3...")
  P3(dt, path)
  print("Part 4...")
  P4(dt, path)
}

P5(dt, "mortality")


get_number_series <- function(path, forecast_horizon = 12) {
  dt <- readRDS(sprintf("%s/ets/eval_%s.rds", path, forecast_horizon))
  dt$dtb %>%
    filter(!startsWith(cluster, "permute")) %>%
    filter(cluster != "average") %>%
    select(cluster, S, batch) %>%
    mutate(cluster = stringi::stri_replace_all_fixed(cluster, "-dr", "")) %>%
    rowwise() %>%
    mutate(S = NROW(S)) %>%
    group_by(cluster) %>%
    summarise(S = mean(S, na.rm = TRUE)) %>%
    mutate(path = path)
}

rbind(get_number_series("mortality"), get_number_series("tourism")) %>%
  tidyr::pivot_wider(values_from = "S", names_from = "path") %>%
  write.csv("manuscript/figures/n_series.csv")



# features
features.compute <- function(data, frequency = frequency) {
  feature_lst <- c(
    "acf_features", "arch_stat", "autocorr_features", "crossing_points", "dist_features",
    "entropy", "heterogeneity", "hurst", "lumpiness", "stability", "pacf_features", "stl_features",
    "unitroot_kpss", "unitroot_pp", "nonlinearity", "max_level_shift", "max_var_shift", "max_kl_shift",
    "holt_parameters", "hw_parameters", "flat_spots"
  )

  ts_features <- tsfeatures::tsfeatures(ts(data$bts, frequency = frequency), features = feature_lst)
  ts_features
}

mortality_feature <- features.compute(readRDS("mortality/ets/batch_144.rds"), frequency = 12)
tourism_features <- features.compute(readRDS("tourism/ets/batch_120.rds"), frequency = 12)

feat_lst <- c("seas_acf1", "hw_parameters_gamma", "trend")

tourism_features %>%
  select(all_of(feat_lst)) %>%
  tidyr::pivot_longer(names_to = "features", cols = everything()) %>%
  mutate(path = "tourism") %>%
  rbind(mortality_feature %>%
    select(all_of(feat_lst)) %>%
    tidyr::pivot_longer(names_to = "features", cols = everything()) %>%
    mutate(path = "mortality")) %>%
  group_by(features, path) %>%
  summarise(value = mean(value)) %>%
  tidyr::pivot_wider(names_from = "path", values_from = "value") %>%
  write.csv("manuscript/figures/features.csv")
