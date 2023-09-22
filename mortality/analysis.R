library(dplyr)
library(tidyr)
# ts and error representator
source("R/construct_hierarchy.R", chdir = TRUE)

rf_method <- "mint"
accuracy_method <- "rmse"


# read all the results
store_all <- NULL
for (batch in 0:11) {
  store_all <- readRDS(paste0("mortality/store_", batch, ".rds"))$output %>%
    output_pre() %>% mutate(batch = batch) %>%
    filter(rf_method == .env$rf_method, accuracy_method == .env$accuracy_method) %>%
    select(-rf_method, -accuracy_method) %>%
    rbind(store_all)
}

store_all <- store_all %>% rowwise() %>%
  mutate(values = list(c(total, bottom))) %>% ungroup() %>%
  select(-total, -bottom)

## remove single random
store_all <- store_all %>% filter(!startsWith(cluster, "random-single"))


# store_all %>% group_by(representator, distance, cluster) %>% tally() %>%
#   pull(n) %>% table()

## combine the results of 12 batches for one clustering method

store_all <- store_all %>% select(-other) %>% nest(values = c(batch, values)) %>%
  mutate_at("values", purrr::map, function(x){
    x <- arrange(x, batch)
    do.call(c, x$values)
  })

## convert the tibble to a data.frame used for comparing ranks
store_mat <- do.call(cbind, store_all$values)
store_rank <- apply(store_mat, 1, rank) %>% t() %>% colMeans()
store_all$rank <- store_rank

## MCB test
mcb <- function(x) {
  tsutils::nemenyi(x, plottype = "vmcb")
}
method_desc <- function(row){
  print(sprintf("representator: %s, distance: %s, cluster: %s", 
                store_all$representator[row],
                store_all$distance[row],
                store_all$cluster[row]))
}

mcb(store_mat)




# representator rank
representator_summarise <- store_all %>% filter(distance != "") %>%
  nest(values = c(representator, values, rank)) %>%
  mutate_at("values", purrr::map, function(x){
    x$rank <- do.call(cbind, x$values) %>% apply(1, rank) %>%
      t() %>% colMeans()
    x %>% select(-values) %>%
      pivot_wider(names_from = "representator", values_from = "rank")
  }) %>%
  unnest("values")


representators <- c("error", "ts", "error.features", "ts.features", "forecast")

## in most cases, error is the best
representator_summarise %>% rowwise() %>%
  mutate(min = representators[which.min(c(error, ts, error.features, ts.features, forecast))]) %>%
  pull(min) %>%
  table()

## only in one case, error is the worst
representator_summarise %>% rowwise() %>%
  mutate(max = representators[which.max(c(error, ts, error.features, ts.features, forecast))]) %>%
  pull(max) %>%
  table()
mcb(representator_summarise %>% select(representators))






# no hierarchy is not good: 391 / 428
rank(store_all$rank)[1]

# natural hierarchy is not bad: 72/428
rank(store_all$rank)[2]

# random versus random natural versus natural

tmp <- store_all %>% filter(distance == "")
tmp$rank <- do.call(cbind, tmp$values) %>% apply(1, rank) %>%
  t() %>% colMeans()

## 1. all the averaged random hierarchies beat no hierarchy
## 2. deep and wide hierarchy wins
## 3. average of more hierarchies does not obtain better results. 
##   Some really bad hierarchies destroy the results. 
##.  Deeper hierarchies have the ability to mitigate the affect of bad series


# cluster rank
cluster_summarise <- store_all %>% filter(distance != "") %>%
  nest(values = c(cluster, values, rank)) %>%
  mutate_at("values", purrr::map, function(x){
    x$rank <- do.call(cbind, x$values) %>% apply(1, rank) %>%
      t() %>% colMeans() %>% rank()
    x %>% select(-values) %>%
      pivot_wider(names_from = "cluster", values_from = "rank")
  }) %>%
  unnest("values")


cluster_summarise %>% select(-c(representator, distance)) %>%
  mcb()

## 1. kmedoids unnested is really bad
## 2. hcluster with single linkage is bad
## 3. nested kmedoids performs really well, which forms a natural-like unbalanced hierarchy.



# distance rank
distance_summarise <- store_all %>% filter(distance != "") %>%
  nest(values = c(distance, values, rank)) %>%
  mutate_at("values", purrr::map, function(x){
    x$rank <- do.call(cbind, x$values) %>% apply(1, rank) %>%
      t() %>% colMeans() %>% rank()
    x %>% select(-values) %>%
      pivot_wider(names_from = "distance", values_from = "rank")
  }) %>%
  unnest("values")
distance_summarise %>% select(-c(representator, cluster)) %>%
  mcb()

distances <- c("euclidean", "dtw", "negcor", "cor", "uncorrelation")
distance_summarise %>% rowwise() %>%
  mutate(max = distances[which.max(c(euclidean, dtw, negcor, cor, uncorrelation))]) %>%
  pull(max) %>%
  table()
distance_summarise %>% rowwise() %>%
  mutate(min = distances[which.min(c(euclidean, dtw, negcor, cor, uncorrelation))]) %>%
  pull(min) %>%
  table()
















