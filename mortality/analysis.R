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
method_desc()





