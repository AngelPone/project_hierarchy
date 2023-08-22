library(dplyr)
# ts and error representator
source("R/construct_hierarchy.R", chdir = TRUE)

store_all <- NULL
for (batch in 0:11) {
  store_all <- readRDS(paste0("mortality/store_", batch, ".rds"))$output %>%
    output_pre() %>% mutate(batch = batch) %>%
    rowwise() %>%
    mutate(bottom_mean = mean(bottom)) %>%
    ungroup() %>%
    rbind(store_all)
}


store_mint_mae <- store_all %>% filter(rf_method == "mint", accuracy_method == "mae") %>%
  select(-rf_method, -accuracy_method, -bottom, -other)


store_mint_rmse <- store_all %>% filter(rf_method == "mint", accuracy_method == "rmse") %>%
  select(-rf_method, -accuracy_method, -bottom, -other)



store_base <- store_all %>% filter(rf_method == "base", accuracy_method == "mae") %>%
  rename(total_base = total, bottom_base = bottom) %>%
  select(batch, total_base, bottom_base)

mae_ratio <- store_all %>% filter(accuracy_method == "mae") %>%
  left_join(store_base, by = "batch") %>% 
  rowwise() %>%
  mutate(total = total/total_base, bottom = list(c(bottom/bottom_base))) %>%
  mutate(bottom_mean = mean(bottom), bottom_max = max(bottom),  bottom_min = min(bottom)) %>%
  ungroup()

store_base_rmse <- store_all %>% filter(rf_method == "base", accuracy_method == "rmse") %>%
  rename(total_base = total, bottom_base = bottom) %>%
  select(batch, total_base, bottom_base)


rmse_ratio <- store_all %>% filter(accuracy_method == "rmse") %>%
  left_join(store_base_rmse, by = "batch") %>% 
  rowwise() %>%
  mutate(total = total/total_base, bottom = list(c(bottom/bottom_base))) %>%
  mutate(bottom_mean = mean(bottom), bottom_max = max(bottom),  bottom_min = min(bottom)) %>%
  ungroup()



# Conclusion 1
mae_ratio %>% group_by(rf_method) %>% summarise(total= mean(total, na.rm=TRUE), bottom_max=mean(bottom_max, na.rm=TRUE), 
                                                bottom_min=mean(bottom_min, na.rm=TRUE), bottom_mean=mean(bottom_mean, na.rm=TRUE))
rmse_ratio %>% group_by(rf_method) %>% summarise(total= mean(total, na.rm=TRUE), bottom_max=mean(bottom_max, na.rm=TRUE), 
                                                bottom_min=mean(bottom_min, na.rm=TRUE), bottom_mean=mean(bottom_mean, na.rm=TRUE))

# the ability to improve all the time series
mae_ratio %>% filter(bottom_max <= 1.05, rf_method != "base") %>% select(representator, distance, cluster, total, bottom_mean, bottom_max, bottom_min, batch) %>%
  arrange(total)

