library(Matrix)

create_new_output <- function() {
  output <- vector("list", 5)
  names(output) <- c("representator", "distance", "cluster", "accuracy", "other")
  output
}

saveResult <- function() {
  saveRDS(list(bfstore = BASEFORECAST_STORE,
               data = data,
               output = output,
               features = FEATURES),
          store_path)
}

add_result <- function(output, representotar, distance, 
                       cluster, accs, other = NULL) {
  if (is.null(output$representator)) {
    output$representator <- representotar
    output$distance <- distance
    output$cluster <- cluster
    output$accuracy <- list(accs)
    output$other <- list(other)
  } else {
    output$representator <- c(output$representator, representotar)
    output$distance <- c(output$distance, distance)
    output$cluster <- c(output$cluster, cluster)
    output$accuracy <- append(output$accuracy, list(accs))
    output$other <- append(output$other, list(other))
  }
  output
}

output_pre <- function(output) {
  as_tibble(output) %>% 
    tidyr::unnest_longer(accuracy, values_to = "accuracy", indices_to = "accuracy_method") %>% 
    tidyr::unnest_longer(accuracy, values_to = "accuracy", indices_to = "rf_method") %>% 
    tidyr::unnest_wider(accuracy)
}

metrics <- c("mae", "rmse")
if (exists("path")) {
  store_path <- sprintf("%s/store_%s.rds", path, batch)
}



