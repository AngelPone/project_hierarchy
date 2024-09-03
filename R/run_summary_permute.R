args <- commandArgs(trailingOnly = TRUE)
path <- args[[1]]

source("R/utils.R")


compute_comb3 <- function(batch) {
  nls <- Filter(function(x) {
    (x$representor != "") & (startsWith(x$cluster, "permute"))
  }, batch$nl)
  combined_rf <- list()
  for (permute in 1:100) {
    nls_ <- Filter(\(x) endsWith(x$cluster, paste0("-", permute)), nls)
    stopifnot(length(nls_) == 12)
    combined_rf[[permute]] <- do.call(abind::abind, list(purrr::map(nls_, function(x){x$rf}), along = 0))
    combined_rf[[permute]] <- apply(combined_rf[[permute]], c(2, 3), mean)
  }
  combined_rf
}

acc <- function() {
  batch_length <- list(mortality=144, tourism=120)
  batch_length <- batch_length[[path]]
  dt <- readRDS(sprintf("%s/data.rds", path))
  dtb <- NULL
  S <- rbind(rep(1, NCOL(dt$S)), diag(NCOL(dt$S)))

  combination1 <- map(readRDS("mortality/combination.rds"), \(x) x[[1]])
  for (batch in 0:batch_length) {
    data <- readRDS(sprintf("%s/batch_%s.rds", path, batch))
    tts <- data$tts %*% t(S)
    bts <- data$bts %*% t(S)
    dtb_ <- nl2tibble(data$nl)
    if (path == "mortality") {
      permute_combination <- compute_comb3(data)
      dtb_ <- add_row(dtb_, representor="", distance="",
                      cluster="combination1", rf=list(combination1[[batch+1]]),
                      S=NULL, other=NULL)
      for (i in 1:100) {
        dtb_ <- add_row(dtb_, representor="", distance="",
                        cluster=paste0("permute-combination1-", i),
                        rf=list(permute_combination[[i]]), S=NULL, other=NULL)
      }
    }
    dtb_[["rmsse"]] <- map_dbl(dtb_$rf, function(g) {
      mean(sapply(1:NCOL(g), \(x) metric.rmsse(tts[,x], g[,x], bts[,x])))
    })
    dtb_$batch <- batch
    dtb <- rbind(dtb, dtb_)
  }
  dtb
}

dt <- acc()

pdf(sprintf("manuscript/figures/%s/natural_vs_pn.pdf", path), width = 8, height = 6)
par(mar=c(4,14,3,2))
natural_hierarchy <- mcb_hierarchy_rmsse(
  dt %>% filter(cluster == "natural"),
  dt %>% filter(startsWith(cluster, "permute-natural")),
  "Natural"
)
dev.off()
rank <- which(names(natural_hierarchy$means) == "Natural")
itl1 <- natural_hierarchy$means - natural_hierarchy$cd / 2
itl2 <- natural_hierarchy$means + natural_hierarchy$cd / 2
sig_better <-
  length(which(itl1 > itl2["Natural"]))
sig_worse <-
  length(which(itl2 < itl1["Natural"]))

write(sprintf("Natural ranks %s in its 100 twins, significantly better than %s, significantly worse than %s",
              rank, sig_better, sig_worse),
      sprintf("manuscript/figures/%s/natural_vs_pn.txt", path))


pdf(sprintf("manuscript/figures/%s/cluster_vs_pc.pdf", path), 8, 6)
par(mar=c(4,14,3,2))
best_ <- readRDS(sprintf("%s/eval_cluster.rds", path))$best_
best_name <- method_name(best_$representor, best_$distance, best_$cluster)
cluster_hierarchy_rmsse <- mcb_hierarchy_rmsse(
  dt %>% filter(
    representor == best_$representor,
    cluster == best_$cluster,
    distance == best_$distance
  ),
  dt %>% filter(
    representor == best_$representor,
    startsWith(cluster, paste0("permute-", best_$cluster)),
    distance == best_$distance
  ),
  method_name(best_$representor, best_$distance, best_$cluster)
)
dev.off()

rank <- which(names(cluster_hierarchy_rmsse$means) == best_name)
itl1 <- cluster_hierarchy_rmsse$means - cluster_hierarchy_rmsse$cd / 2
itl2 <- cluster_hierarchy_rmsse$means + cluster_hierarchy_rmsse$cd / 2
sig_better <- length(which(itl1 > itl2[best_name]))
sig_worse <- length(which(itl2 < itl1[best_name]))

write(sprintf("The best cluster ranks %s in its 100 twins, significantly better than %s, significantly worse than %s",
              rank, sig_better, sig_worse),
      sprintf("manuscript/figures/%s/cluster_vs_pc.txt", path))



if (path == "mortality") {
  # calculate combination of twin hierarchies
  pdf(sprintf("manuscript/figures/%s/comb_vs_pc.pdf", path), 8, 6)
  par(mar=c(4,14,3,2))
  mcb <- mcb_hierarchy_rmsse(
    dt %>% filter(cluster == "combination1"),
    dt %>% filter(startsWith(cluster, "permute-combination1")),
    "Combination"
  )
  dev.off()
}

