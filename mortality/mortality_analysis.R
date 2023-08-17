library(dplyr)
# ts and error representator
output <- readRDS("mortality/ts_error.rds")
d1 <- output_pre(output) %>% rowwise() %>%
  mutate(bottom = mean(bottom), num_aggregate_series = as.character(NROW(other$S))) %>%
  select(-other) %>% 
  arrange(total, bottom)

output <- readRDS("mortality/output_base.rds")
d2 <- output_pre(output) %>% 
  rowwise() %>%
  mutate(bottom = mean(bottom), 
         num_aggregate_series = ifelse(is.null(other), "0", paste0(other$n_clusters, " * 20 "))) %>%
  select(-other)


d <- d1 %>% add_row(d2) %>%
  ungroup() %>%
  arrange(total)

