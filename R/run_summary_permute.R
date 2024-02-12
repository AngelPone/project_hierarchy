args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

path <- "tourism"

dt <- readRDS(sprintf("%s/%s/eval.rds", path, "ets"))
dt_orig <- readRDS(sprintf("%s/%s/batch_0.rds", path, "ets"))
dt_names <- readRDS(sprintf("%s/names.rds", path))
source("R/run_nls.R")
library(dplyr)
library(ggplot2)
library(Matrix)
library(tsutils)


df_bench <- dt$dtb %>%
  filter(!startsWith(cluster, "permute"))

df_bench %>% rowwise() %>% mutate(rmsse=mean(rmsse)) %>%
  group_by(representor, distance, cluster) %>%
  summarise(rmsse=mean(rmsse)) %>%
  arrange(rmsse)

df_bench %>% rowwise() %>% mutate(rmsse=mean(rmsse)) %>%
  select(representor, distance, cluster, batch, rmsse) %>%
  mutate(method = paste0(representor, "-", distance, "-", cluster)) %>%
  select(method, batch, rmsse) %>%
  tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
  select(-batch) %>%
  nemenyi(plottype = "vmcb")



clusters <- dt$dtb %>% select(cluster, distance, representor) %>%
  filter(!startsWith(cluster, "permute-")) %>%
  filter(distance != "") %>%
  unique()


# Natural

representor <- ""
distance <- ""
cluster <- "natural"
permutesd <- dt$dtb %>% filter(representor==.env$representor, distance==.env$distance,
                               startsWith(cluster, paste0("permute-", .env$cluster)))
bench <- dt$dtb %>% filter(representor==.env$representor, distance==.env$distance,
                           cluster==.env$cluster)
all_rmsse <- NULL
for (batch in unique(dt$dtb$batch)) {
  bench_batch <- bench %>% filter(batch == .env$batch)
  permutesd_batch <- permutesd %>% filter(batch == .env$batch) %>% 
    arrange(cluster)
  permutesd_rmsse <- permutesd_batch %>% rowwise() %>%
    mutate(rmsse = mean(rmsse)) %>% pull(rmsse)
  bench_rmsse <- mean(bench_batch$rmsse[[1]])
  all_rmsse <- rbind(all_rmsse, c(bench_rmsse, permutesd_rmsse))
}
colnames(all_rmsse) <- c("Natural", 1:100)
test_ <- tsutils::nemenyi(all_rmsse, plottype = "vmcb")
print(paste("rank of natural in the sorted random hierarchies: ", which(colnames(test_$intervals) == "Natural")))
test_$means

# cluster

all_rmsse2 <- NULL
bench_rmsse <- clusters %>% left_join(dt$dtb) %>%
  rowwise() %>% mutate(rmsse=mean(rmsse)) %>%
  select(representor, distance, cluster, batch, rmsse)

best_ <- bench_rmsse %>% group_by(representor, distance, cluster) %>% summarise(rmsse=mean(rmsse)) %>%
  arrange(rmsse)
best_ <- best_[1,] %>% select(representor, distance, cluster)
for (batch in unique(bench_rmsse$batch)) {
  best_batch <- best_ %>% left_join(bench_rmsse, by = c("representor", "distance", "cluster")) %>%
    filter(batch == .env$batch)
  permute_best_ <- dt$dtb %>% filter(
    representor == best_$representor[[1]],
    distance == best_$distance[[1]],
    startsWith(cluster, paste0("permute-", best_$cluster[[1]])),
    batch == .env$batch
  ) %>% rowwise() %>% mutate(rmsse = mean(rmsse)) %>%
    pull(rmsse)
  all_rmsse2 <- rbind(all_rmsse2, c(best_batch$rmsse[[1]], permute_best_))
  stopifnot(NROW(permute_best_) == 100)
}
colnames(all_rmsse2) <- c("bestcluster", 1:100)
test_ <- tsutils::nemenyi(all_rmsse2, plottype = "vmcb")
print(paste("rank of best cluster method in its sorted random hierarchies: ", which(colnames(test_$intervals) == "bestcluster")))
test_$means


if (path == "mortality") {
  
  
  
  # averaging
  all_rmsse3 <- NULL
  for (batch in unique(bench_rmsse$batch)) {
    store_path <- sprintf("%s/%s/batch_%s.rds", path, bfmethod, batch)
    data <- readRDS(store_path)
    data_tibble <- nl2tibble(data$nl) %>% 
      filter(!startsWith(cluster, "random"),
             representor != "",
             !startsWith(cluster, "hcluster-random"))
    
    data_tibble$permute <- sapply(data_tibble$cluster, function(x){
      if (!startsWith(x, "permute")) {
        return (0)
      }
      x <- strsplit(x, "-")[[1]]
      x <- x[length(x)]
      x
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
        sapply(x$rf, function(g) { g[["mint"]] }, simplify = "array") %>%
          apply(c(1,2), mean)
      })
    
    tts <- cbind(rowSums(data$tts), data$tts)
    bts <- cbind(rowSums(data$bts), data$bts)
    
    rmsses <- sapply(iterators::iter(avg$rf), function(f) {
      sapply(1:NCOL(f), function(x){
        metric.rmsse(tts[,x], f[,x], bts[,x])
      }) %>% mean()
    })
    all_rmsse3 <- rbind(all_rmsse3, c(rmsses[1], sort(rmsses[2:101])))
  }
  colnames(all_rmsse3) <- c("average", 1:100)
  test_ <- tsutils::nemenyi(all_rmsse3, plottype = "none")
  print(paste("rank of average cluster method in its sorted random hierarchies: ", which(colnames(test_$intervals) == "average")))
  
  
  # nemenyi(cbind(all_rmsse[,1], all_rmsse2[,1], all_rmsse3[,1]), plottype = "vmcb")
  # 
  # mat <- cbind(all_rmsse3, all_rmsse2[,1], all_rmsse[,1])
  # colnames(mat) <- c("average", 1:100, "bestcluster", "bench")
  # test_ <- nemenyi(mat, plottype = "vmcb")
  # print(paste("rank of average cluster method in its sorted random hierarchies: ", which(colnames(test_$intervals) == "average")))
  # print(paste("rank of best cluster method in average random hierarchies: ", which(colnames(test_$intervals) == "bestcluster")))
  # print(paste("rank of natural method in average random hierarchies: ", which(colnames(test_$intervals) == "bench")))
  # 
  
}

