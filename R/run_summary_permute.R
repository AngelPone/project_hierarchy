source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)
library(tsutils)

forecast_horizon <- 1
source("R/metrics.R")
source("R/expr_utils.R")

method_name <- function(representor, distance, cluster) {
  representor <- strsplit(representor, "-")[[1]][1]
  representor <- ifelse(!is.na(representor), switch(representor, error = "ER-", error.features = "ERF-",
                                                    ts = "TS-", ts.features = "TSF-"), "")
  distance <- ifelse(distance == "","",
                     switch(distance, euclidean="-EUC", dtw="-DTW")
  )
  cluster <- strsplit(cluster, "-")[[1]][1]
  cluster <- ifelse(
    !is.na(cluster),
    switch(cluster, natural="Natural", Kmedoids="ME", hcluster="HC", base="Base", "average"="Average"),
    "Two-level"
  )
  paste0(representor, cluster, distance)
}




mcb_series_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>% arrange(batch) %>%
    pull(rmsse) %>%
    do.call(c, .)
  rand_rmsse <- 
    rand %>% arrange(batch, cluster) %>%
    select(cluster, rmsse, batch) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    mutate_at("rmsse", purrr::map, function(g){
      g <- g %>% arrange(batch)
      list(do.call(c, g$rmsse))
    }) %>%
    tidyr::unnest(rmsse) %>%
    pull(rmsse) %>%
    do.call(cbind, .)
  all_rmsse <- cbind(orig_rmsse, rand_rmsse)
  colnames(all_rmsse) <- c(name, 1:100)
  nemenyi(all_rmsse, plottype = "vmcb", target = name)
}

mcb_hierarchy_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>% rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch) %>%
    pull(rmsse)
  
  rand_rmsse <- rand %>% rowwise() %>% 
    mutate(rmsse = mean(rmsse)) %>%
    arrange(batch, cluster) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    pull(rmsse) %>%
    lapply(function(x) { x %>% arrange(batch) %>% pull(rmsse) }) %>%
    do.call(cbind, .)
  
  all_rmsse <- cbind(orig_rmsse, rand_rmsse)
  colnames(all_rmsse) <- c(name, 1:100)
  nemenyi(all_rmsse, plottype = "vmcb", target = name)
}


rmsse_benchmarks <- NULL
rank_natural_tbl <- list()
rank_cluster_tbl <- list()


series_mcb <- function(dt) {
  a <- dt %>% select(method, batch, rmsse) %>%
    tidyr::nest(rmsse = -"method") %>%
    mutate_at("rmsse", purrr::map, function(x){
      x %>% arrange(batch) %>%
        pull(rmsse) %>%
        do.call(c, .)
    })
  
  dat <- do.call(cbind, a$rmsse)
  colnames(dat) <- a$method
  dat
}

# natural vs two-level
P1 <- function(dt, path){
  a <- dt$dtb %>% filter(cluster %in% c("natural", "")) %>%
    mutate(method = ifelse(cluster=="natural", "Natural", "Two-level")) %>%
    select(method, rmsse, batch)
  
  b <- a %>% rowwise() %>% mutate(rmsse = mean(rmsse))
  
  pdf(sprintf("manuscript/figures/hierarchy_rmsse/%s/P1_natural_vs_twolevel_h%s.pdf", path, forecast_horizon))
  b %>% tidyr::pivot_wider(names_from = "method", values_from = "rmsse") %>%
    select(`Two-level`, Natural) %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()
  
  pdf(sprintf("manuscript/figures/series_rmsse/%s/P1_natural_vs_twolevel_h%s.pdf", path, forecast_horizon))
  a %>% series_mcb() %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()
  
  P1_table[[path]] <<- b %>% group_by(method) %>%
    summarise(rmsse = mean(rmsse))
  
  if (length(P1_table) == 2) {
    P1_table <- P1_table[["tourism"]] %>%
      mutate(dataset = "tourism") %>%
      rbind(P1_table[["mortality"]] %>% mutate(dataset="mortality")) %>%
      tidyr::pivot_wider(names_from = "dataset", values_from = "rmsse")
    write.csv(P1_table, sprintf("manuscript/figures/hierarchy_rmsse/P1_rmsse_h%s.csv", forecast_horizon))
  }
}

# natural vs its randomization
P2 <- function(dt, path) {
  # MCB Test based on single series rmsse
  pdf(sprintf("manuscript/figures/series_rmsse/%s/P2_natural_vs_pn_h%s.pdf", path, forecast_horizon),
      width = 8, height = 6)
  natural_series <- mcb_series_rmsse(
    dt$dtb %>% filter(cluster == "natural"),
    dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
    "Natural"
  )
  dev.off()
  
  # MCB Test based on hierarchy rmsse
  pdf(sprintf("manuscript/figures/hierarchy_rmsse/%s/P2_natural_vs_pn_h%s.pdf", path, forecast_horizon),
      width = 8, height = 6)
  natural_hierarchy <- mcb_hierarchy_rmsse(
    dt$dtb %>% filter(cluster == "natural"),
    dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
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
  rmsse_tbl <- P1_table[[path]]
  
  P2_table[[path]] <<- c(rank, rmsse_tbl$rmsse[which(rmsse_tbl$method == "Natural")], 
                         sig_better, sig_worse)
  if (length(P2_table) == 2) {
    data.frame(P2_table) %>%
      write.csv(sprintf("manuscript/figures/hierarchy_rmsse/P2_rmsse_h%s.csv", forecast_horizon))
  }
}


P3 <- function(dt, path) {
  bench_rmsse <- dt$dtb %>%
    filter(!startsWith(cluster, "permute"), cluster != "average") %>%
    rowwise() %>%
    mutate(rmsse=mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster="base") %>%
            rowwise() %>% mutate(rmsse=mean(rmsse)))
  
  pdf(sprintf("manuscript/figures/hierarchy_rmsse/%s/P3_mcb_benchmarks_h%s.pdf", path, forecast_horizon), 
      width = 8, height = 6)
  par(mex=1.1)
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
    summarise(rmsse = mean(rmsse) * 100, .groups = "drop") %>%
    filter(representor != "") %>%
    arrange(rmsse)
  
  best_ <- best_[1,]
  
  best_name <- method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])
  
  pdf(sprintf("manuscript/figures/series_rmsse/%s/P3_cluster_vs_pc_h%s.pdf", path, forecast_horizon),
      width = 8, height = 6)
  cluster_series_rmsse <- mcb_series_rmsse(
    dt$dtb %>% filter(representor == best_$representor[[1]],
                      cluster == best_$cluster[[1]],
                      distance == best_$distance[[1]]),
    dt$dtb %>% filter(representor == best_$representor[[1]],
                      startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
                      distance == best_$distance[[1]]),
    best_name
  )
  dev.off()
  
  
  pdf(sprintf("manuscript/figures/hierarchy_rmsse/%s/P3_cluster_vs_pc_h%s.pdf", path, forecast_horizon),
      width = 8, height = 6)
  cluster_hierarchy_rmsse <- mcb_hierarchy_rmsse(
    dt$dtb %>% filter(representor == best_$representor[[1]],
                      cluster == best_$cluster[[1]],
                      distance == best_$distance[[1]]),
    dt$dtb %>% filter(representor == best_$representor[[1]],
                      startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
                      distance == best_$distance[[1]]),
    method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])
  )
  dev.off()
  
  
  # calculate rmsse
  bench_rmsse <- bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    group_by(method) %>%
    summarise(rmsse = mean(rmsse) * 100) %>%
    mutate(rmsse = round(rmsse, digits = 3))
  
  itl1 <- cluster_hierarchy_rmsse$means - cluster_series_rmsse$cd / 2
  itl2 <- cluster_hierarchy_rmsse$means + cluster_series_rmsse$cd / 2
  sig_better <- 
    which(itl1 > itl2[best_name])
  sig_worse <-
    which(itl2 < itl1[best_name])
  P3_table[[path]] <<- bench_rmsse %>% mutate(dataset=path)
  P3_rank[[path]] <<- c(which(names(cluster_hierarchy_rmsse$means) == best_name),
                        P3_table[[path]]$rmsse[which(P3_table[[path]]$method == best_name)],
                        length(sig_better),
                        length(sig_worse))
  
  if (length(P3_table) == 2) {
    methods <- c("Base", "Two-level","Natural", "TS-HC-EUC", "TS-HC-DTW", "TS-ME-EUC", "TS-ME-DTW",
                 "TSF-HC-EUC", "TSF-ME-EUC",
                 "ER-HC-EUC", "ER-HC-DTW", "ER-ME-EUC", "ER-ME-DTW",
                 "ERF-HC-EUC", "ERF-ME-EUC")
    a <- do.call(rbind, P3_table) %>%
      tidyr::pivot_wider(names_from = "dataset", values_from = "rmsse")
    a[match(methods, a$method),] %>%
      write.csv(sprintf("manuscript/figures/hierarchy_rmsse/P3_rmsse_h%s.csv", forecast_horizon))
    
    data.frame(P3_rank) %>%
      write.csv(sprintf("manuscript/figures/hierarchy_rmsse/P3_rank_h%s.csv", forecast_horizon))
  }
}

P4 <- function(dt, path) {
  
  bench_rmsse <- dt$dtb %>%
    filter(!startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(rmsse=mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster="base") %>%
            rowwise() %>% mutate(rmsse=mean(rmsse))) %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse)
  
  P4_table[[path]] <<- bench_rmsse %>% group_by(method) %>%
    summarise(rmsse = round(mean(rmsse) * 100, digits = 3))
  
  if (path == "mortality") {
    pdf(sprintf("manuscript/figures/hierarchy_rmsse/%s/P4_average_vs_pa_h%s.pdf", path, forecast_horizon), height = 6, width = 8)
      test_ <- mcb_hierarchy_rmsse(
        dt$dtb %>% filter(representor == "", cluster == "average",distance == ""),
        dt$dtb %>% filter(representor=="", distance=="", startsWith(cluster, "permute-average")),
        "Average"
      )
    dev.off()
    
    itl1 <- test_$means - test_$cd / 2
    itl2 <- test_$means + test_$cd / 2
    sig_better <- 
      which(itl1 > itl2["Average"])
    sig_worse <- 
      which(itl2 < itl1["Average"])
    rmsse_ <- P4_table[[path]]$rmsse[which(P4_table[[path]]$method == "Average")]
    P4_rank[[path]] <-
      c(rmsse_, which(names(test_$means) == "Average"), length(sig_better), length(sig_worse))
    data.frame(P4_rank) %>% write.csv(sprintf("manuscript/figures/hierarchy_rmsse/P4_rank_h%s.csv", forecast_horizon))
    
    pdf(sprintf("manuscript/figures/series_rmsse/%s/P4_average_vs_pa_h%s.pdf", path, forecast_horizon), height = 6, width = 8)
    test_ <- mcb_series_rmsse(
      dt$dtb %>% filter(representor == "", cluster == "average",distance == ""),
      dt$dtb %>% filter(representor=="", distance=="", startsWith(cluster, "permute-average")),
      "Average"
    )
    dev.off()
  }
  
  pdf(sprintf("manuscript/figures/hierarchy_rmsse/%s/P4_benchmarks_h%s.pdf", path, forecast_horizon), height = 6, width = 8)
  bench_rmsse %>%  
    tidyr::pivot_wider(names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()
  

  methods <- c("Base", "Two-level","Natural", "TS-HC-EUC", "TS-HC-DTW", "TS-ME-EUC", "TS-ME-DTW",
               "TSF-HC-EUC", "TSF-ME-EUC",
               "ER-HC-EUC", "ER-HC-DTW", "ER-ME-EUC", "ER-ME-DTW",
               "ERF-HC-EUC", "ERF-ME-EUC", "Average")
  P4_table[[path]] <<-
    P4_table[[path]][match(methods, P4_table[[path]]$method),] %>%
    mutate(dataset = path)
  
  if (length(P4_table) == 2) {
    rbind(P4_table[[1]], P4_table[[2]]) %>%
      tidyr::pivot_wider(names_from = "dataset", values_from = "rmsse") %>%
      write.csv(sprintf("manuscript/figures/hierarchy_rmsse/P4_rmsse_h%s.csv", forecast_horizon))
  }
}

for (forecast_horizon in c(1, 12)) {
  P1_table <- list()
  P2_table <- list()
  P3_table <- list()
  P3_rank <- list()
  P4_table <- NULL
  P4_rank <- list()
  for (path in c("tourism", "mortality")) {
      dt <- readRDS(sprintf("%s/ets/eval_%s.rds", path, forecast_horizon))
      print("Part 1...")
      P1(dt, path)
      print("Part 2...")
      P2(dt, path)
      print("Part 3...")
      P3(dt, path)
      print("Part 4...")
      P4(dt, path)
  }
}
