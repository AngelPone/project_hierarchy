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
    switch(cluster, natural="Natural", Kmedoids="ME", hcluster="HC", base="Base"),
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

for (path in c("tourism", "mortality")) {
  dt <- readRDS(sprintf("%s/%s/eval.rds", path, "ets"))
  # MCB Test based on single series rmsse
  # pdf(sprintf("manuscript/figures/%s_natural_mcb_single.pdf", path),width = 8, height = 6)
  # natural_series <- mcb_series_rmsse(
  #   dt$dtb %>% filter(cluster == "natural"),
  #   dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
  #   "Natural"
  # )
  # dev.off()
  # which(names(natural_series$means) == "Natural")
  
  # MCB Test based on hierarchy rmsse
  pdf(sprintf("manuscript/figures/%s_natural_mcb_whole.pdf", path),width = 8, height = 6)
  natural_hierarchy <- mcb_hierarchy_rmsse(
    dt$dtb %>% filter(cluster == "natural"),
    dt$dtb %>% filter(startsWith(cluster, "permute-natural")),
    "Natural"
  )
  dev.off()
  
  
  # base, natural, cluster rmsse
  bench_rmsse <- dt$dtb %>%
    filter(!startsWith(cluster, "permute")) %>%
    rowwise() %>%
    mutate(rmsse=mean(rmsse)) %>%
    select(representor, distance, cluster, batch, rmsse) %>%
    rbind(dt$base %>% mutate(representor = "", distance = "", cluster="base") %>%
            rowwise() %>% mutate(rmsse=mean(rmsse)))
  
  
  pdf(sprintf("manuscript/figures/%s_benchmark_mcb_whole.pdf", path),width = 8, height = 6, family = "Helvetica")
  par(mex=1.1)
  bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
    select(-batch) %>%
    tsutils::nemenyi(plottype = "vmcb")
  dev.off()
  
  
  clusters <- dt$dtb %>% select(cluster, distance, representor) %>%
    filter(!startsWith(cluster, "permute-")) %>%
    filter(distance != "") %>%
    unique()
  
  # cluster
  
  best_ <- bench_rmsse %>%
    group_by(representor, distance, cluster) %>%
    summarise(rmsse = mean(rmsse) * 100, .groups = "drop") %>%
    filter(representor != "") %>%
    arrange(rmsse)
  
  best_ <- best_[1,]
  
  best_name <- method_name(best_$representor[[1]], best_$distance[[1]], best_$cluster[[1]])
  
  # pdf(sprintf("manuscript/figures/%s_cluster_mcb_single.pdf", path),width = 8, height = 6)
  # cluster_series_rmsse <- mcb_series_rmsse(
  #   dt$dtb %>% filter(representor == best_$representor[[1]],
  #                     cluster == best_$cluster[[1]],
  #                     distance == best_$distance[[1]]),
  #   dt$dtb %>% filter(representor == best_$representor[[1]],
  #                     startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
  #                     distance == best_$distance[[1]]),
  #   best_name
  # )
  # dev.off()
  # cluster_series_rmsse$means
  
  
  pdf(sprintf("manuscript/figures/%s_cluster_mcb_whole.pdf", path),width = 8, height = 6)
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
  
  bench_rmsse <- bench_rmsse %>%
    rowwise() %>%
    mutate(method = method_name(representor, distance, cluster)) %>%
    select(method, batch, rmsse) %>%
    group_by(method) %>%
    summarise(rmsse = mean(rmsse) * 100) %>%
    mutate(rmsse = round(rmsse, digits = 3)) %>%
    arrange(rmsse)
  
  if (path == "tourism") {
    rmsse_benchmarks <- bench_rmsse %>% rename(tourism = "rmsse")
  } else {
    rmsse_benchmarks <- bench_rmsse %>% 
      rename(mortality = "rmsse") %>%
      left_join(rmsse_benchmarks, by = "method")
  }
  rank_natural_tbl[[path]] <- c(
    bench_rmsse$rmsse[which(bench_rmsse$method == "Natural")],
    which(names(natural_hierarchy$means) == "Natural"))
  rank_cluster_tbl[[path]] <- c(
    bench_rmsse$rmsse[which(bench_rmsse$method == best_name)],
    which(names(cluster_hierarchy_rmsse$means) == best_name))
  
  if (path == "mortality") {
    # averaging
    rank_average_tbl <- list()
    all_rmsse3 <- NULL
    for (batch in unique(dt$dtb$batch)) {
      store_path <- sprintf("%s/%s/batch_%s.rds", path, "ets", batch)
      data <- readRDS(store_path)
      data_tibble <- nl2tibble(data$nl) %>% 
        filter(representor != "")
      
      data_tibble$permute <- sapply(data_tibble$cluster, function(x){
        if (!startsWith(x, "permute")) {
          return (0)
        }
        x <- strsplit(x, "-")[[1]]
        x <- x[length(x)]
        as.integer(x)
      })
      data_tibble$cluster <- sapply(data_tibble$cluster, function(x){
        if (!startsWith(x, "permute")) {
          return (x)
        }
        split_x <- strsplit(x, "-")[[1]]
        x <- stringi::stri_replace_all_fixed(x, "permute-", "")
        x <- stringi::stri_replace_all_fixed(x, paste0("-", split_x[length(split_x)]), "")
        x
      })
      
      avg <- data_tibble %>% select(representor, cluster, distance, rf, permute) %>%
        arrange(permute) %>%
        group_by(permute) %>%
        tidyr::nest(rf = -"permute") %>%
        mutate_at("rf", purrr::map, function(x) {
          sapply(x$rf, function(g) { g[["mint"]][1:forecast_horizon,,drop=FALSE] }, simplify = "array") %>%
            apply(c(1,2), mean)
        })
      if (is.null(dim(data$tts))){
        tts <- c(sum(data$tts), data$tts)
        tts <- matrix(tts, nrow=1)
      } else { 
        tts <- cbind(rowSums(data$tts), data$tts)[1:forecast_horizon,,drop=FALSE]
      }
      bts <- cbind(rowSums(data$bts), data$bts)
      
      rmsses <- sapply(iterators::iter(avg$rf), function(f) {
        sapply(1:NCOL(f), function(x){
          metric.rmsse(tts[,x], f[,x], bts[,x])
        }) %>% mean()
      })
      all_rmsse3 <- rbind(all_rmsse3, rmsses)
    }
    colnames(all_rmsse3) <- c("Average", 1:100)
    
    pdf("manuscript/figures/mortality_average_mcb_whole.pdf", height = 6, width = 8)
    test_ <- nemenyi(all_rmsse3, plottype = "vmcb", target = "Average")
    dev.off()
    rank_average_tbl[["mortality"]] <-
      c(round(colMeans(all_rmsse3)[1]*100, digits=3), which(names(test_$means) == "Average"))
  }
}


write.csv(rmsse_benchmarks, "manuscript/figures/cluster_rmsse.csv")

rank_natural_tbl <- data.frame(rank_natural_tbl)
row.names(rank_natural_tbl) <- c("RMSSE", "Rank")
write.csv(rank_natural_tbl, "manuscript/figures/rank_natural.csv")

rank_cluster_tbl <- data.frame(rank_cluster_tbl)
row.names(rank_cluster_tbl) <- c("RMSSE", "Rank")
write.csv(rank_cluster_tbl, "manuscript/figures/rank_cluster.csv")


rank_average_tbl <- data.frame(rank_average_tbl)
row.names(rank_average_tbl) <- c("RMSSE", "Rank")
write.csv(rank_average_tbl, "manuscript/figures/rank_average.csv")
