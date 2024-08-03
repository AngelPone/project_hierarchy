library(Matrix, quietly = TRUE)

create_new_output <- function() {
  output <- vector("list", 5)
  names(output) <- c("representator", "distance", "cluster", "accuracy", "other")
  output
}



output_pre <- function(output) {
  as_tibble(output) %>% 
    tidyr::unnest_longer(accuracy, values_to = "accuracy", indices_to = "accuracy_method") %>% 
    tidyr::unnest_longer(accuracy, values_to = "accuracy", indices_to = "rf_method") %>% 
    tidyr::unnest_wider(accuracy)
}



mcb_ <- function(orig, type, plot_type, target,
                 file_path, mar = NULL, conf.level) {
  orig_ <- orig %>%
    select(method, rmsse, batch) %>%
    tidyr::nest(rmsse = -c("method")) %>%
    mutate_at("rmsse", purrr::map, function(g) {
      g <- g %>% arrange(batch)
      if (type == "hierarchy") {
        return(list(sapply(g$rmsse, mean)))
      } else {
        return(list(do.call(c, g$rmsse)))
      }
    }) %>%
    tidyr::unnest(rmsse)
  
  rmsse <- do.call(cbind, orig_$rmsse)
  colnames(rmsse) <- orig_$method
  pdf(file_path, width=8, height = 6)
  if (plot_type == "custom") {
    stopifnot(!is.null(mar))
    par(mar = mar)
    nemenyi(rmsse, plottype = "vmcb", target = target, conf.level = conf.level)
  }
  if (plot_type == "orig") {
    tsutils::nemenyi(rmsse, plottype = "vmcb", conf.level = conf.level)
  }
  dev.off()
}
mcb_series_rmsse <- function(orig, rand, name) {
  orig_rmsse <- orig %>%
    arrange(batch) %>%
    pull(rmsse) %>%
    do.call(c, .)
  rand_rmsse <-
    rand %>%
    arrange(batch, cluster) %>%
    select(cluster, rmsse, batch) %>%
    tidyr::nest(rmsse = -"cluster") %>%
    mutate_at("rmsse", purrr::map, function(g) {
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
