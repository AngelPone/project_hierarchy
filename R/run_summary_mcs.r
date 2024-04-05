source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)
library(tsutils)
library(MCS)

source("R/metrics.R")
source("R/expr_utils.R")


cl <- parallel::makeCluster(8)

method_name <- function(representor, distance, cluster) {
  representor <- strsplit(representor, "-")[[1]][1]
  representor <- ifelse(!is.na(representor), switch(representor,
    error = "ER-",
    error.features = "ERF-",
    ts = "TS-",
    ts.features = "TSF-"
  ), "")
  distance <- ifelse(distance == "", "",
    switch(distance,
      euclidean = "EUC",
      dtw = "DTW"
    )
  )
  cluster <- strsplit(cluster, "-")[[1]][1]
  cluster <- ifelse(
    !is.na(cluster),
    switch(cluster,
      natural = "Natural",
      Kmedoids = "-ME",
      hcluster = "-HC",
      base = "Base",
      "average" = "Combination"
    ),
    "Two-level"
  )
  paste0(representor, distance, cluster)
}

mcs_hierarchy_rmsse <- function(orig, rand, name, alpha) {
  orig_rmsse <- orig %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch) %>%
    pull(rmsse)

  rand_rmsse <- rand %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch, cluster) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    pull(rmsse) %>%
    lapply(function(x) {
      x %>%
        arrange(batch) %>%
        pull(rmsse)
    }) %>%
    do.call(cbind, .)

  all_rmsse <- cbind(orig_rmsse, rand_rmsse)
  colnames(all_rmsse) <- c(name, 1:100)
  MCSprocedure(all_rmsse, alpha = alpha, cl=cl)
}


# natural vs two-level
P1 <- function(dt, path, alpha) {
  a <- dt$dtb %>%
    filter(cluster %in% c("natural", "")) %>%
    mutate(method = ifelse(cluster == "natural", "Natural", "Two-level")) %>%
    select(rmsse, batch, method) %>%
    rbind(dt$base %>% mutate(method = 'Base'))

  b <- a %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse))

  b %>%
    tidyr::pivot_wider(names_from = "method", values_from = "rmsse") %>%
    select(`Two-level`, Natural, Base) -> b
  
  MCSprocedure(b, alpha = alpha, cl = cl)
}

# natural vs its randomization
P2 <- function(dt, path, alpha) {
  mcs_hierarchy_rmsse(
    dt$dtb %>% filter(cluster == "natural"),
    dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
    "Natural",
    alpha
  )
}


P3 <- function(dt, path, alpha) {
  bench_rmsse <- dt$dtb %>%
    filter(!startsWith(cluster, "permute"), cluster != "average") %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster = "base") %>%
      rowwise() %>% mutate(rmsse = mean(rmsse)))

  bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    MCSprocedure(alpha = alpha, cl=cl) -> benchmark

  best_ <- bench_rmsse %>%
    group_by(representor, distance, cluster) %>%
    summarise(rmsse = mean(rmsse), .groups = "drop") %>%
    filter(representor != "") %>%
    arrange(rmsse)

  best_ <- best_[1, ]

  best_name <- method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])

  mcs <- mcs_hierarchy_rmsse(
    dt$dtb %>% filter(
      representor == best_$representor[[1]],
      cluster == best_$cluster[[1]],
      distance == best_$distance[[1]]
    ),
    dt$dtb %>% filter(
      representor == best_$representor[[1]],
      startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
      distance == best_$distance[[1]]
    ),
    method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]]),
    alpha
  )
  
  return(list(benchmark=benchmark, permutation=mcs))
}

P4 <- function(dt, path, alpha) {
  bench_rmsse <- dt$dtb %>%
    filter(!startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster = "base") %>%
      rowwise() %>% mutate(rmsse = mean(rmsse))) %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse)
  
  bench_rmsse %>%
    tidyr::pivot_wider(names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    MCSprocedure(alpha = alpha, cl=cl) -> benchmark

  test_ <- NULL
  if (path == "mortality") {
    
    test_ <- mcs_hierarchy_rmsse(
      dt$dtb %>% filter(representor == "", cluster == "average", distance == ""),
      dt$dtb %>% filter(representor == "", distance == "", startsWith(cluster, "permute-average")),
      "Combination",
      alpha
    )
  }
  
  return(list(benchmark=benchmark, permutation=test_))
}


output85 <- list()
for (path in c("tourism", "mortality")) {
  dt <- readRDS(sprintf("%s/ets/eval_%s.rds", path, forecast_horizon))
  print("Part 1...")
  output85[[path]] <- list()
  output85[[path]]$set1 <- P1(dt, path, 0.15)
  print("Part 2...")
  output85[[path]]$set2 <- P2(dt, path, 0.15)
  print("Part 3...")
  output85[[path]]$set3 <- P3(dt, path, 0.15)
  print("Part 4...")
  output85[[path]]$set4 <- P4(dt, path, 0.15)
}

get_mcs_info <- function(res, target) {
  c(length(res@Info$model.names), target %in% res@Info$model.names)
}




best_ <- list(tourism = "TSF-EUC-HC", mortality = "TS-DTW-HC")
summ <- list()
for (path in c("mortality", "tourism")) {
  best_clustering <-
  summ[[path]] <- rbind(get_mcs_info(output[[path]]$set1, "Natural"),
                        get_mcs_info(output[[path]]$set2, "Natural"),
                        get_mcs_info(output[[path]]$set3$benchmark, "Natural"),
                        get_mcs_info(output[[path]]$set3$permutation, best_[[path]]),
                        get_mcs_info(output[[path]]$set4$benchmark, "Combination"),
                        get_mcs_info(output[[path]]$set4$benchmark, "Natural"))
  rownames(summ[[path]]) <- c("P1", "Natural vs PN", "Cluster vs Natural", "Cluster vs PC", "Combination", "Combination - Natural")
  if (path == "mortality") {
    summ[[path]] <- rbind(summ[[path]], get_mcs_info(output[[path]]$set4$permutation, "Combination"))
    rownames(summ[[path]]) <- c("P1", "Natural vs PN", "Cluster vs Natural", "Cluster vs PC", "Combination", "Combination - Natural", "Combination vs PC")
  }
}
summ








