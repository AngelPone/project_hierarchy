get_number_series <- function(path, forecast_horizon = 12) {
  dt <- readRDS(sprintf("%s/eval_cluster.rds", path, forecast_horizon))
  dt$dtb %>%
    filter(!startsWith(cluster, "permute")) %>%
    filter(!cluster %in% c("", "combination1", "combination2", "base")) %>%
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

mortality_feature <- features.compute(readRDS("mortality/batch_144.rds"), frequency = 12)
tourism_features <- features.compute(readRDS("tourism/batch_120.rds"), frequency = 12)

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
