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

rank_compare <- function(dt, methods_df = NULL, measure="rmse", rf_method = "mint", filter_random = FALSE) {
  if (is.null(methods_df)) {
    if (filter_random) {
      methods_random <- dt$dtb %>% select(representor, distance, cluster) %>% 
        filter(cluster %in% c("random-15-50", "hcluster-random-50", "random-natural")) %>% 
        unique()
    } else {
      methods_random <- dt$dtb %>% select(representor, distance, cluster) %>% 
        filter(grepl(cluster, "random", fixed = TRUE)) %>% 
        unique()
    }
    methods_dr <- dt$dtb %>% select(representor, distance, cluster) %>% filter(endsWith(cluster, "-dr")) %>% unique()
    methods_natural <- dt$dtb %>% select(representor, distance, cluster) %>%
      filter(cluster %in% c("base", "", "natural")) %>% unique()
    methods_df <- Kmedoids_clusterN(dt) %>% filter(`0` != length(dt$base$rmse) | is.na(`0`)) %>% 
      select(representor, distance) %>% mutate(cluster = "Kmedoids") %>%
      rbind(methods_dr, methods_random, methods_natural)
  }
  dtb <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    select(representor, distance, cluster, all_of(measure), batch) %>%
    tidyr::nest(values = c("batch", measure)) %>%
    mutate_at("values", purrr::map, function(g) {
      g <- arrange(g, "batch") %>% pull(measure)
      do.call(c, lapply(g, function(x) x[[rf_method]]))
    })
  avg_n <- dt$dtb %>% right_join(methods_df, by = c("representor", "distance", "cluster")) %>%
    rowwise() %>%
    mutate(n = ifelse(is.null(S), 0, NROW(S))) %>%
    group_by(representor, distance, cluster) %>% summarise(n = mean(n), .groups = "drop")
  metric_mat <- do.call(cbind, dtb$values)
  dt_base <- do.call(c, dt$base[[measure]])
  metric_mat <- cbind(metric_mat, dt_base)
  dtb <- dtb %>% mutate(method_name = paste0(cluster,
                                    ifelse(representor == "", "", 
                                           paste0("-", representor, "-", distance)))) %>%
    mutate(method_name = ifelse(method_name == "", "no-cluster", method_name))

  colnames(metric_mat) <- c(dtb$method_name, "base")
  mcbtest <- tsutils::nemenyi(metric_mat, plot = "vmcb")
  lowest_rank <- min(mcbtest$means)
  mcbmeans <- mcbtest$means[c(dtb$method_name, "base")]
  dtb %>% select(representor, distance, cluster) %>% 
    add_row(representor = "", distance = "", cluster = "base") %>%
    mutate(avgrank = mcbmeans) %>%
    mutate(sigworse = (mcbmeans > (lowest_rank + mcbtest$cd))) %>%
    right_join(avg_n, by = c("representor", "cluster", "distance")) %>%
    arrange(avgrank)
}


Kmedoids_clusterN <- function(dt, cluster = "Kmedoids", mode="S") {
  dt$dtb %>% filter(cluster == .env$cluster) %>% 
    rowwise() %>%
    mutate(n_cluster = ifelse(is.null(S), 0, ifelse(mode == "S", NROW(S)))) %>%
    select(representor, distance, batch, n_cluster) %>% ungroup() %>%
    group_by(representor, distance, n_cluster) %>%
    summarise(count=n(), .groups = "drop") %>%
    arrange(n_cluster) %>%
    pivot_wider(id_cols = c("representor", "distance"), names_from = "n_cluster", values_from = "count")
}

Kmedoids_gap <- function(dt) {
  dt$dtb %>% filter(cluster == "Kmedoids-dr") %>%
    rowwise() %>%
    mutate(n_cluster = ifelse(is.null(S), 0, NROW(S))) %>%
    mutate(n_Kmedoids = other$n_cluster) %>%
    mutate(n_gap = other$gap) %>%
    ungroup() %>%
    select(-c(S, other, rmsse, mase, rmse, mae))
}

inspect_silhouette <- function(dt_orig, representor) {
  source("R/representator.R")
  distance_mat <- dt_orig$distance[[representor]][["euclidean"]]
  nl <- nl2tibble(dt_orig$nl)
  idx <- which(nl$representor == representor & nl$cluster == "Kmedoids-dr")
  
  print(summary(silhouette(pam(distance_mat, k=nl$other[[idx]]$n_cluster,diss=TRUE))))
  if (nl$other[[idx]]$gap > 1) {
    print(summary(silhouette(pam(distance_mat, k=nl$other[[idx]]$gap,diss=TRUE))))
  }
}

visualizeDistance <- function(dt, representor, distance, cluster = "Kmedoids") {
  distance_mat <- dt$distance[[representor]][[distance]]
  d2scaled <- cmdscale(distance_mat, 2, eig = TRUE)
  d2scaled <- tibble(x=d2scaled$points[,1], y=d2scaled$points[,2])
  tmpdt <- nl2tibble(dt$nl) %>% filter(representor == .env$representor, 
                                   distance == .env$distance,
                                   cluster == .env$cluster)
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
  gap <- tmpdt$other[[1]]$gap
  p + ggtitle(sprintf("%s %s %s %s", representor, distance, 
                      paste0("nKmedoids = ", tmpdt$other[[1]]$n_cluster),
                      ifelse(is.null(gap), "", paste0("gap = ", gap))))
}
  



visualizeGroup <- function(dt, representor, distance, type = c("bts", "resid"), cluster = "Kmedoids", names = NULL) {
  tmpdt <- nl2tibble(dt$nl) %>% filter(cluster == .env$cluster, representor == .env$representor,
                                       distance == .env$distance)
  cluster_S <- tmpdt$S[[1]]
  
  distance <- dt$distance[[representor]][[distance]]
  n_cluster <- tmpdt$other[[1]]$n_cluster
  pam_x <- pam(distance, k=n_cluster, diss = TRUE)
  sil <- summary(silhouette(pam_x))$clus.avg.width
  sil <- sil[which(sil > 0.05)]
  
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
    
    plot(series, main = sprintf("%s Group = %s, group silhouette = %s", type, i, sil[i]))
  }
}
