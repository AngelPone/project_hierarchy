
metric.rmse <- function(obs, pred, hist) {
  sqrt(mean((obs - pred)^2))
}

metric.mae <- function(obs, pred, hist) {
  mean(abs(obs - pred))
}

metric.rmsse <- function(obs, pred, hist) {
  sqrt(mean((obs - pred)^2) / mean(diff(hist, 12)^2))
}

metric.mase <- function(obs, pred, hist) {
  mean(abs(obs - pred)) / mean(abs(diff(hist, 12)))
}


# evaluate.hts <- function(x, metrics, type = c("base", "rf", "nl", "average")) {
#   type <- match.arg(type)
#   tts <- x$tts %*% t(x$S)
#   allts <- x$bts %*% t(x$S)
#   
#   output <- list()
#   for (metric in metrics) {
#     metric_func <- get(paste0("metric.", metric))
#     if (type == "base") {
#       accs <- sapply(1:NCOL(tts), function(i){
#         metric_func(tts[, i], x$basef[, i], allts[,i])
#       })
#       accs <- list(base = list(total = accs[1], bottom = accs[2:length(accs)]))
#     } else if (type == "rf") {
#       accs <- lapply(x$rf, function(rf) {
#         tmp <- sapply(1:NCOL(tts), function(i){
#           metric_func(tts[, i], rf[, i], allts[,i])
#         })
#         tmp <- list(total = tmp[1], bottom = tmp[2:length(tmp)])
#         tmp
#       })
#     } else if(type == "nl") {
#       accs <- list()
#       for (i in seq_along(x$nl)) {
#         accs[[i]] <- lapply(x$nl[[i]]$rf, function(rf){
#           tmp <- sapply(1:NCOL(tts), function(i){
#             metric_func(tts[, i], rf[, i], allts[,i])
#           })
#           tmp <- list(total = tmp[1], bottom = tmp[2:length(tmp)])
#           tmp
#         })
#       }
#       if (length(accs) == 1){
#         accs <- accs[[1]]
#       }
#     } else if(type == "average") {
#       accs <- list()
#       for (rf_method in names(x$nl[[1]]$rf)) {
#         rf <- x$nl[[1]]$rf[[rf_method]]
#         for (i in 2:length(x$nl)) {
#           rf <- rf + x$nl[[i]]$rf[[rf_method]]
#         }
#         rf <- rf/length(x$nl)
#         tmp_accs <- sapply(1:NCOL(tts), function(i){
#           metric_func(tts[, i], rf[, i], allts[,i])
#         })
#         accs[[rf_method]] <- list(total = tmp_accs[1], bottom = tmp_accs[2:length(tmp_accs)])
#       }
#     }
#     output[[metric]] <- accs
#   }
#   output
# }
