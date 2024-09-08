args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]

source("R/utils.R")




hts.eval <- function(df, tts, bts, S, combination, base) {
  tts <- tts %*% t(S)
  bts <- bts %*% t(S)
  
  df <- add_row(df, representor="", cluster="combination1", 
                distance="", S=NULL, 
                rf=list(combination[[1]]), other=NULL) 
  
  df <- add_row(df, representor = "", cluster="combination2", S=NULL, 
                distance = "",
                rf=list(combination[[2]]), other=NULL)
  
  df[["rmsse"]] <- lapply(iterators::iter(df$rf), function(g) {
    c <- g[, 2:(m+1),drop=FALSE]
    c <- c %*% t(S)
    sapply(1:NCOL(c), function(x) { metric.rmsse(tts[,x], c[,x], bts[,x]) })
  })
  
  
  basef_middle <- df$other[df$cluster == "natural"][[1]]$basef
  basef <- cbind(base[,1], basef_middle, base[,2:(m+1)])
  basef_rmsse <- sapply(1:NCOL(bts),
                        function(x) {
                          metric.rmsse(tts[,x], basef[,x], bts[,x])
                        } )
  
  df <- add_row(df, representor="", cluster="base", distance="",
                S=NULL, rf=NULL, rmsse=list(basef_rmsse))
  df
}

dtb <- NULL

combination <- readRDS(sprintf("%s/combination.rds", path))
for (batch in 0:(batch_length-1)) {
  store_path <- sprintf("%s/batch_%s.rds", path, batch)
  data <- readRDS(store_path)
  data_tibble <- nl2tibble(data$nl) %>%
    filter(!startsWith(cluster, "permute"))
  data_tibble <- hts.eval(data_tibble, data$tts, data$bts, dt$S, combination[[batch+1]], data$basef)
  
  data_tibble <- data_tibble %>% select(-rf, -other) %>% mutate(batch = batch)
  dtb <- rbind(dtb, data_tibble)
}


#' mcb all cluster hierarchies, grouped, base, natural, two-level
bench_rmsse <- dtb %>%
  filter(!startsWith(cluster, "permute"), cluster != "combination1") %>%
  rowwise() %>%
  mutate(rmsse = mean(rmsse[c(1, (n-m+1):n)])) %>%
  select(representor, distance, cluster, batch, rmsse) %>%
  ungroup()

stopifnot(NROW(bench_rmsse) == batch_length * 16)

jpeg(sprintf("figures/Figure4_%s_mcb_benchmarks.jpg", path), height = 600/1.3, width = 800/1.3)
par(mex = 1.1)
bench_rmsse %>%
  rowwise() %>%
  mutate(method = method_name(representor, distance, cluster)) %>%
  select(method, batch, rmsse) %>%
  tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
  select(-batch) %>%
  tsutils::nemenyi(plottype = "vmcb")
dev.off()


#' mcb all cluster hierarchies, grouped, base, natural, two-level and combination
bench_rmsse <- dtb %>%
  filter(!startsWith(cluster, "permute")) %>%
  rowwise() %>%
  mutate(rmsse = mean(rmsse[c(1, (n-m+1):n)])) %>%
  select(representor, distance, cluster, batch, rmsse) %>%
  ungroup()

stopifnot(NROW(bench_rmsse) == batch_length * 17)

jpeg(sprintf("figures/Figure12_%s_mcb_combination.jpg", path), height = 600/1.3, width = 800/1.3)
par(mex = 1.1)
bench_rmsse %>%
  rowwise() %>%
  mutate(method = method_name(representor, distance, cluster)) %>%
  select(method, batch, rmsse) %>%
  tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "rmsse") %>%
  select(-batch) %>%
  tsutils::nemenyi(plottype = "vmcb")
dev.off()


methods <- c(
  "Base", "Two-level", "Natural", "Grouped",
  "TS-EUC-ME", "ER-EUC-ME", "TSF-EUC-ME", "ERF-EUC-ME",
  "TS-EUC-HC", "ER-EUC-HC", "TSF-EUC-HC", "ERF-EUC-HC",
  "TS-DTW-ME", "TS-DTW-HC", "ER-DTW-ME", "ER-DTW-HC",
  "Combination"
)

#' save table to csv 
#' table only for total and bottom level
csv_ <- bench_rmsse %>% rowwise() %>%
  mutate(method = method_name(representor, distance, cluster)) %>%
  ungroup() %>%
  group_by(method, representor, distance, cluster) %>% summarise(rmsse = mean(rmsse))
csv_[match(methods, csv_$method), c("method", "rmsse")] %>%
  write.csv(sprintf("figures/Table3_%s_Table8.csv", path))

csv_ %>% filter(representor != "") %>%
  arrange(rmsse) -> csv_
best_cluster_method <- list(representor=csv_$representor[1],
                            cluster=csv_$cluster[1],
                            distance=csv_$distance[1])

#' table for total, middle and bottom level, used in supplymentary materials
bench_rmsse <- dtb %>%
  filter(!startsWith(cluster, "permute")) %>%
  rowwise() %>%
  mutate(top = rmsse[1], middle=mean(rmsse[2:(n-m)]), bottom=mean(rmsse[(n-m+1):n]),
         method = method_name(representor, distance, cluster)) %>%
  ungroup() %>%
  select(method, batch, top, middle, bottom)

jpeg(sprintf("figures/FigureS1_%s_mcb_combination_middle.jpg", path), height = 600/1.3, width = 800/1.3)
par(mex = 1.1)
bench_rmsse %>%
  select(method, batch, middle) %>%
  tidyr::pivot_wider(id_cols = "batch", names_from = "method", values_from = "middle") %>%
  select(-batch) %>%
  tsutils::nemenyi(plottype = "vmcb")
dev.off()


bench_rmsse <- bench_rmsse %>%
  group_by(method) %>%
  summarise(across(c(top, middle, bottom), \(x) mean(x)))

bench_rmsse[match(methods, bench_rmsse$method),] %>%
  write.csv(sprintf("figuresTableS1_%s_rmsse_threelevel.csv", path))

saveRDS(list(dtb=dtb,
             best_ = best_cluster_method), sprintf("%s/eval_cluster.rds", path))

