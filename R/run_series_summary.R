source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)
library(tsutils)

source("R/metrics.R")
source("R/expr_utils.R")

# natural vs two-level
P1 <- function(dt, path, type = "series", conf.level = 0.95) {
  
  
  path <- sprintf("manuscript/figures/%s_rmsse/%s/P1_natural_vs_twolevel_h%s.pdf", type, path, forecast_horizon)
  
  mcb_(dt %>% filter(method %in% c("Natural", "Two-level", "Base")), type = type, 
       plot_type = "orig", file_path = path, conf.level = conf.level)
}

# natural vs its randomization
P2 <- function(dt, path, type = "series", conf.level = 0.95) {
  path <- sprintf("manuscript/figures/%s_rmsse/%s/P2_natural_vs_pn_h%s.pdf", type, path, forecast_horizon)
  
  mcb_(dt %>% 
         filter(cluster == "natural" | startsWith(cluster, "permute-natural")),
       type = type, plot_type = "custom", file_path = path, mar=c(4,14,3,2), 
       target = "Natural", conf.level = conf.level)
}



P3 <- function(dt, path, type = "series", conf.level = 0.95) {
  
  file_path <- sprintf("manuscript/figures/%s_rmsse/%s/P3_mcb_benchmarks_h%s.pdf", type, path, forecast_horizon)
  mcb_(dt %>% 
         filter(!startsWith(cluster, "permute"), cluster != "average"),
       type = type, plot_type = "orig", file_path = file_path, conf.level = conf.level)
  
  
  bench_rmsse <- dt %>%
    filter(representor != "", !startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    select(representor, distance, cluster, method, batch, rmsse) %>%
    group_by(representor, distance, cluster, method) %>%
    summarise(rmsse = mean(rmsse), .groups = "drop") %>%
    arrange(rmsse)

  best_ <- bench_rmsse[1, ]
  file_path <- sprintf("manuscript/figures/%s_rmsse/%s/P3_cluster_vs_pc_h%s.pdf", type, path, forecast_horizon)
  mcb_(dt %>% 
         filter(method == best_$method[1] | startsWith(cluster, paste0("permute-", best_$cluster[1]))) %>%
         filter(representor == best_$representor[1]),
       type = type, plot_type = "custom", file_path = file_path, mar=c(4,14,3,2), target = best_$method[1], conf.level = conf.level)
}

P4 <- function(dt, path, type = "series", conf.level = 0.95) {
  if (path == "mortality") {
    file_path <- sprintf("manuscript/figures/%s_rmsse/%s/P4_average_vs_pa_h%s.pdf", type, path, forecast_horizon)
    mcb_(dt %>% filter(cluster == "average" | startsWith(cluster, "permute-average")),
         type, "custom", "Combination", file_path, mar=c(4,14,3,2), conf.level = conf.level)
  }
  
  file_path <- sprintf("manuscript/figures/%s_rmsse/%s/P4_benchmarks_h%s.pdf", type, path, forecast_horizon)
  mcb_(dt %>% filter(!startsWith(cluster, "permute")),
       type = type, "orig", "Combination", file_path, conf.level = conf.level)
}

forecast_horizon <- 12
type <- "hierarchy"
conf.level <- 0.9

for (path in c("tourism", "mortality")) {
  dt <- readRDS(sprintf("%s/ets/eval_%s.rds", path, forecast_horizon))
  dt <- dt$dtb %>% rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster = "Base", 
                             method = "Base", S = list(NULL), other = list(NULL))) %>%
    ungroup()
  print("Part 1...")
  P1(dt, path, type, conf.level)
  print("Part 2...")
  P2(dt, path, type, conf.level)
  print("Part 3...")
  P3(dt, path, type, conf.level)
  print("Part 4...")
  P4(dt, path, type, conf.level)
}

