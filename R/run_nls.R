# Calculate average RMSE to compare with mint paper
nl2tibble <- function(x) {
  output <- vector("list", 6)
  names(output) <- c("representor", "cluster", "distance", "S", "rf", "other")
  
  for (i in seq_along(x)) {
    for (n in names(x[[i]])) {
      if (is.character(x[[i]][[n]])) {
        output[[n]] <- append(output[[n]], x[[i]][[n]])
      } else {
        output[[n]] <- append(output[[n]], list(x[[i]][[n]]))
      }
    }
  }
  tibble::as_tibble(output)
}

avg_measure_fn <- function(dt, rf_method = "mint", metric = "rmse") {
  n <- length(dt$base$rmse[[1]])
  base <- dt$base %>% pull(metric) %>% 
    do.call(rbind, .) %>%
    colMeans()
  dt$dtb %>% 
    nest(!!metric := -c("representor", "distance", "cluster")) %>% 
    mutate_at(metric, purrr::map, function(x){
      x <- arrange(x, batch)
      x %>% pull(metric) %>%
        lapply(function (z) z[[rf_method]]) %>%
        do.call(rbind, .) %>% colMeans()
    }) %>%
    add_row(representor = "", distance = "", cluster="base", !!metric  := list(base)) %>%
    rowwise() %>%
    mutate(total = get(metric)[1], bottom = mean(get(metric)[2:n])) %>%
    arrange(total)
}

rank_compare <- function(dt, methods_df = NULL, measure="rmse", rf_method = "mint") {
  if (is.null(methods_df)) {
    methods_df <- Kmedoids_clusterN(dt) %>% filter(`0` != length(dt$base$rmse) | is.na(`0`)) %>% 
      select(representor, distance) %>% mutate(cluster = "Kmedoids") %>%
      rbind(dt$dtb %>% select(representor, distance, cluster) %>% filter(distance == "") %>% unique())
  }
  dtb <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    select(representor, distance, cluster, all_of(measure), batch) %>%
    tidyr::nest(values = c("batch", measure)) %>%
    mutate_at("values", purrr::map, function(g) {
      g <- arrange(g, "batch") %>% pull(measure)
      do.call(c, lapply(g, function(x) x[[rf_method]]))
    })
  metric_mat <- do.call(cbind, dtb$values)
  dt_base <- do.call(c, dt$base[[measure]])
  metric_mat <- cbind(metric_mat, dt_base)
  dtb <- dtb %>% mutate(method_name = paste0(cluster,
                                    ifelse(representor == "", "", 
                                           paste0("-", representor, "-", distance)))) %>%
    mutate(method_name = ifelse(method_name == "", "no-cluster", method_name))

  colnames(metric_mat) <- c(dtb$method_name, "base")
  mcbtest <- tsutils::nemenyi(metric_mat, plot = "vline")
  lowest_rank <- min(mcbtest$means)
  mcbmeans <- mcbtest$means[c(dtb$method_name, "base")]
  dtb %>% select(representor, distance, cluster) %>% 
    add_row(representor = "", distance = "", cluster = "base") %>%
    mutate(avgrank = mcbmeans) %>%
    mutate(sigworse = (mcbmeans > (lowest_rank + mcbtest$cd))) %>%
    arrange(avgrank)
}

# Compare rank
# rank_compare <- function(dt, measure, rf_method, col, col1, col2) {
#   col1_choice <- dt %>% filter(representor != "") %>% 
#     pull({{col1}}) %>% unique()
#   col2_choice <- dt %>% filter(representor != "") %>% 
#     pull({{col2}}) %>% unique()
#   
#   output <- vector("list", 4)
#   names(output) <- c("best", "sigbest", "worst", "sigworse")
#   
#   ch1 <- rep(col1_choice, length(col2_choice))
#   ch2 <- rep(col2_choice, each = length(col1_choice))
#   
#   
#   output <- foreach(choice1 = ch1, choice2 = ch2, .packages = "dplyr") %do% {
#     tmpdt <- dt %>% filter({{col1}} == choice1, {{col2}} == choice2) %>%
#       select(-S) %>%
#       select(batch, representor, cluster, distance, all_of(measure)) %>%
#       tidyr::nest(data = c("batch", measure))
#     
#     # tmpdt <- rbind(tmpdt1, tmpdt2)
#     tmpnames <- c(tmpdt %>% pull({{col}}))
#     
#     dat_mat <- do.call(cbind, lapply(tmpdt$data, function(g) {
#       lapply(g[[measure]], function(e) {
#         if (rf_method %in% names(e)) { return(e[[rf_method]]) }
#         e
#       }) %>% do.call(c, .)
#     }))
#     colnames(dat_mat) <- tmpnames
#     mcb_res <- tsutils::nemenyi(dat_mat, plot = "none")
#     
#     is.bestsignificant <- function(x){
#       x$mean[1] + x$cd/2 < x$mean[2] - x$cd/2
#     }
#     which.best <- function(x){
#       names(x$means)[1]
#     }
#     
#     which.worst <- function(x){
#       names(x$means)[length(x$means)]
#     }
#     
#     is.worstsignificant <- function(x){
#       x$mean[1] + x$cd/2 < x$mean[x$k] - x$cd/2
#     }
#     
#     list(best=which.best(mcb_res), sigbest=is.bestsignificant(mcb_res),
#       worse=which.worst(mcb_res), sigworse=is.worstsignificant(mcb_res))
#   }
#   
#   list2tibble <- function(x, idx) {
#     sapply(x, function(x) {x[[idx]]} )
#   }
#   output <- list(best = list2tibble(output, 1), sigbest = list2tibble(output, 2),
#                  worst = list2tibble(output, 3), sigworse = list2tibble(output, 4))
#   
#   list(best = as_tibble(output) %>% group_by(best) %>% summarise(n(), sigbest = sum(sigbest)),
#        worst = as_tibble(output) %>% group_by(worst) %>% summarise(n(), sigworst = sum(sigworse)))
# }


Kmedoids_clusterN <- function(dt) {
  dt$dtb %>% filter(cluster == "Kmedoids") %>% 
    rowwise() %>%
    mutate(n_cluster = ifelse(is.null(S), 0, NROW(S))) %>%
    select(representor, distance, batch, n_cluster) %>% ungroup() %>%
    group_by(representor, distance, n_cluster) %>%
    summarise(count=n(), .groups = "drop") %>%
    arrange(n_cluster) %>%
    pivot_wider(id_cols = c("representor", "distance"), names_from = "n_cluster", values_from = "count")
}


visualizeDistance <- function(dt, representor, distance) {
  distance_mat <- dt$distance[[representor]][[distance]]
  d2scaled <- cmdscale(distance_mat, 2, eig = TRUE)
  d2scaled <- tibble(x=d2scaled$points[,1], y=d2scaled$points[,2])
  tmpdt <- nl2tibble(dt$nl) %>% filter(representor == .env$representor, 
                                   distance == .env$distance,
                                   cluster == "Kmedoids")
  S <- tmpdt$S[[1]]
  
  p <- ggplot(d2scaled, mapping = aes(x=x,y=y))
  if (!is.null(S)) {
    grpvec <- vector("integer", NCOL(S))
    sapply(1:NROW(S), function(x) { grpvec[which(S[x,] == 1)] <<- x })
    grpvec[tmpdt$other[[1]]$medoids] <- "medoids"
    d2scaled$cluster <-  grpvec
    p <- p + geom_point(aes(color = cluster, group = cluster), data = d2scaled)
  } else {
    p <- p + geom_point()
  }
  p + ggtitle(sprintf("%s %s", representor, distance))
}
  



visualizeGroup <- function(dt, representor, distance, type = c("bts", "resid"), cluster = "Kmedoids", names = NULL) {
  tmpdt <- nl2tibble(dt$nl) %>% filter(cluster == .env$cluster, representor == .env$representor,
                                       distance == .env$distance)
  cluster_S <- tmpdt$S[[1]]
  type <- match.arg(type)
  for (i in 1:NROW(cluster_S)) {
    idx <- which(cluster_S[i,] == 1)[1:min(8, sum(cluster_S[i,]))]
    cnames <- names[idx]
    if (type == "resid") { idx <- idx + 1 }
    series <- dt[[type]][,idx]
    if (!is.null(cnames)) {
      colnames(series) <- cnames
    }
    series <- ts(series, frequency = 12)
    
    plot(series, main = sprintf("%s Group = %s", type, i))
  }
}
