args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]
bfmethod <- args[[2]]

dt <- readRDS(sprintf("%s/%s/eval.rds", path, bfmethod))
dtb <- dt$base
dt <- dt$dtb
library(dplyr)

# Calculate average RMSE to compare with mint paper
avg_measure_fn <- function(dt) {
  output <- list()
  
  for (measure in c("rmse", "mae", "rmsse", "mase")) {
    for (rf_method in c("ols", "wlss", "wlsv", "mint")) {
      tmp <- dt %>% filter(cluster == "natural") %>% pull(measure) %>% 
        lapply(function(x) {x[[rf_method]]} ) %>%
        do.call(rbind, .) %>%
        colMeans()
      output[[rf_method]] <- c(output[[rf_method]], c(tmp[1], mean(tmp[2:length(tmp)])))
    }
    output[["metric"]] <- c(output[["metric"]], rep(measure, 2))
  }
  
  base_output <- list()
  for (measure in c("rmse", "mae", "rmsse", "mase")) {
    tmp <- dtb %>% pull(measure) %>% 
      do.call(rbind, .) %>%
      colMeans()
    base_output[["base"]] <- c(base_output[["base"]], c(tmp[1], mean(tmp[2:length(tmp)])))
    base_output[["metric"]] <- c(base_output[["metric"]], rep(measure, 2))
  }
  list(natural = as_tibble(output), base = as_tibble(base_output))
}

(avg_measure <- avg_measure_fn(dt))


# Compare rank
rank_compare <- function(dt, measure, rf_method, col, col1, col2) {
  col1_choice <- dt %>% filter(representor != "") %>% 
    pull({{col1}}) %>% unique()
  col2_choice <- dt %>% filter(representor != "") %>% 
    pull({{col2}}) %>% unique()
  
  output <- vector("list", 4)
  names(output) <- c("best", "sigbest", "worst", "sigworse")
  
  ch1 <- rep(col1_choice, length(col2_choice))
  ch2 <- rep(col2_choice, each = length(col1_choice))
  
  # tmpdt1 <- rbind(dt %>% filter(cluster == "natural"), 
  #                 dt %>% filter(cluster == ""), 
  #                 dt %>% filter(startsWith(cluster, "random"))) %>%
  #   select(-S) %>%
  #   rbind(dtb %>% mutate(representor = "", cluster = "base", distance = ""), .) %>%
  #   select(batch, representor, cluster, distance, all_of(measure)) %>%
  #   tidyr::nest(data = c("batch", measure))
  # tmpdt1 <- NULL
  
  
  output <- foreach(choice1 = ch1, choice2 = ch2, .packages = "dplyr") %do% {
    tmpdt <- dt %>% filter({{col1}} == choice1, {{col2}} == choice2) %>%
      select(-S) %>%
      select(batch, representor, cluster, distance, all_of(measure)) %>%
      tidyr::nest(data = c("batch", measure))
    
    # tmpdt <- rbind(tmpdt1, tmpdt2)
    tmpnames <- c(tmpdt %>% pull({{col}}))
    
    dat_mat <- do.call(cbind, lapply(tmpdt$data, function(g) {
      lapply(g[[measure]], function(e) {
        if (rf_method %in% names(e)) { return(e[[rf_method]]) }
        e
      }) %>% do.call(c, .)
    }))
    colnames(dat_mat) <- tmpnames
    mcb_res <- tsutils::nemenyi(dat_mat, plot = "none")
    
    is.bestsignificant <- function(x){
      x$mean[1] + x$cd/2 < x$mean[2] - x$cd/2
    }
    which.best <- function(x){
      names(x$means)[1]
    }
    
    which.worst <- function(x){
      names(x$means)[length(x$means)]
    }
    
    is.worstsignificant <- function(x){
      x$mean[1] + x$cd/2 < x$mean[x$k] - x$cd/2
    }
    
    list(best=which.best(mcb_res), sigbest=is.bestsignificant(mcb_res),
      worse=which.worst(mcb_res), sigworse=is.worstsignificant(mcb_res))
  }
  
  list2tibble <- function(x, idx) {
    sapply(x, function(x) {x[[idx]]} )
  }
  output <- list(best = list2tibble(output, 1), sigbest = list2tibble(output, 2),
                 worst = list2tibble(output, 3), sigworse = list2tibble(output, 4))
  
  list(best = as_tibble(output) %>% group_by(best) %>% summarise(n(), sigbest = sum(sigbest)),
       worst = as_tibble(output) %>% group_by(worst) %>% summarise(n(), sigworst = sum(sigworse)))
}

# rank_compare(dt, "mae", "mint", "representor", "cluster", "distance")


  
