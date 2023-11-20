path <- "tourism"
bfmethod <- "ets"
dt <- readRDS(sprintf("%s/%s/eval.rds", path, bfmethod))
dt_orig <- readRDS(sprintf("%s/%s/batch_0.rds", path, bfmethod))
dt_names <- readRDS(sprintf("%s/names.rds", path))
source("R/run_nls.R")
library(dplyr)
library(ggplot2)

rank_compare_summary <- function(dt) {
  
  methods_random <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(cluster %in% c("random-15-10", "random-natural-10")) %>% 
    unique()
  methods_random$name <- ifelse(methods_random$cluster == "random-15-10", "FC-R", "FC-N")
  
  methods_dr <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(endsWith(cluster, "-dr")) %>% unique()
  
  switch_name <- function(x){
    x <- strsplit(x, "-")[[1]][1]
    switch(x, ts="TS", ts.features="TSF", error="ER", error.features = "ERF")
  }
  
  methods_dr$name <- paste0(sapply(methods_dr$representor, switch_name), "-",
                            ifelse(methods_dr$cluster == "hcluster-dr", "HC", "ME"))
  methods_dtw <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(distance == "dtw", representor %in% c("ts", "error")) %>% unique()
  methods_dtw$name <- paste0(ifelse(methods_dtw$representor == "ts", "TS", "ER"), "-",
                             ifelse(methods_dtw$cluster == "hcluster", "HC", "ME"), "-DTW")
  methods_natural <- dt$dtb %>% select(representor, distance, cluster) %>%
    filter(cluster %in% c("base", "", "natural", "cluster-average")) %>% unique()
  methods_natural$name <- ifelse(methods_natural$cluster == "base", "BASE",
                                 ifelse(methods_natural$cluster == "", "two-level",
                                        ifelse(methods_natural$cluster == "natural", "natural", "FC-C")))
  
  methods_df <- rbind(methods_random, methods_dr, methods_dtw, methods_natural)

  dtb <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    select(representor, distance, cluster, rmse, batch, name) %>%
    tidyr::nest(values = c("batch", "rmse")) %>%
    mutate_at("values", purrr::map, function(g) {
      g <- arrange(g, "batch") %>% pull("rmse")
      do.call(c, lapply(g, function(x) x[["mint"]]))
    })

  metric_mat <- do.call(cbind, dtb$values)
  dt_base <- do.call(c, dt$base[["rmse"]])
  metric_mat <- cbind(metric_mat, dt_base)
  
  colnames(metric_mat) <- c(dtb$name, "base")
  tsutils::nemenyi(metric_mat, plot = "vmcb")
  # lowest_rank <- min(mcbtest$means)
  # mcbmeans <- mcbtest$means[c(dtb$method_name, "base")]
  # dtb %>% select(representor, distance, cluster) %>% 
  #   add_row(representor = "", distance = "", cluster = "base") %>%
  #   mutate(avgrank = mcbmeans) %>%
  #   mutate(sigworse = (mcbmeans > (lowest_rank + mcbtest$cd))) %>%
  #   right_join(avg_n, by = c("representor", "cluster", "distance")) %>%
  #   arrange(avgrank)
}


pdf(sprintf("manuscript/figures/%s_mcb.pdf", path), width = 10, height = 6)
par(mar = c(5.1, 10, 4.1, 2.1))
rank_compare_summary(dt)
dev.off()





rank_compare_number <- function(dt) {
  
  methods_random <- dt$dtb %>% select(representor, distance, cluster) %>% 
    filter(cluster %in% paste0(rep(c("random-15-", "random-natural-"), 3), rep(c(10, 20, 50), each=2))) %>% 
    unique()
  name <- ifelse(startsWith(methods_random$cluster, "random-15"), "FC-R", "FC-N")
  name2 <- sapply(methods_random$cluster, function(x){
    paste0("-", strsplit(x, "-")[[1]][3])
  })
  name2[name2=="-10"] <- ""
  methods_random$name <- paste0(name, name2)

  methods_natural <- dt$dtb %>% select(representor, distance, cluster) %>%
    filter(cluster == "cluster-average") %>% unique()
  methods_natural$name <- "FC-C"
  
  methods_df <- rbind(methods_random, methods_natural)
  
  dtb <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    select(representor, distance, cluster, rmse, batch, name) %>%
    tidyr::nest(values = c("batch", "rmse")) %>%
    mutate_at("values", purrr::map, function(g) {
      g <- arrange(g, "batch") %>% pull("rmse")
      do.call(c, lapply(g, function(x) x[["mint"]]))
    })
  
  metric_mat <- do.call(cbind, dtb$values)
  colnames(metric_mat) <- dtb$name
  tsutils::nemenyi(metric_mat, plot = "vmcb")
}

pdf(sprintf("manuscript/figures/%s_number_mcb.pdf", path), width = 6, height = 6)
par(mar = c(5.1, 10, 4.1, 2.1))
rank_compare_number(dt)
dev.off()




