
source("R/utils.R", chdir = T)

compute_comb1 <- function(batch) {
  nls <- Filter(function(x) {
    (x$representor != "") & (!startsWith(x$cluster, "permute"))
  }, batch$nl)
  
  stopifnot(length(nls) == 12)
  combined_rf <- do.call(abind::abind, list(purrr::map(nls, "rf"),
                             along = 0))
  apply(combined_rf, c(2, 3), mean)
}



compute_comb2 <- function(batch) {
  nls <- Filter(function(x) {
    (x$representor != "") & (!startsWith(x$cluster, "permute"))
  }, batch$nl)
  stopifnot(length(nls) == 12)

  S <- do.call(rbind, map(nls, "S"))

  allts <- batch$bts %*% t(as.matrix(S))

  basef <- future_map(as.list(iterators::iter(allts, by="column")),
                      \(x) f.ets(x, 12, 12))
  resids <- do.call(cbind, map(basef, \(x) as.numeric(x$resid)))
  basef <- do.call(cbind, map(basef, \(x) as.numeric(x$basef)))
  basef <- cbind(batch$basef[,1], basef, batch$basef[,2:NCOL(batch$basef)])
  resids <- cbind(batch$resid[,1], resids, batch$resid[,2:NCOL(batch$resid)])

  m <- NCOL(S)
  C <- rep(1, m)
  S <- rbind(C, S)
  rf <- FoReco::csrec(basef, S, comb="shr", res=resids)

  rf[,c(1, (NCOL(rf)-m+1):NCOL(rf))]
}

for (path in c("mortality")) {
  files <- list.files(paste0(path, ""))
  files <- files[startsWith(files, "batch")]
  output <- list()
  for (batch in 1:length(files)) {
    print(sprintf("batch %s ...", batch))
    store_path <- sprintf("%s/batch_%s.rds", path, batch-1)
    dt <- readRDS(store_path)
    comb1 <- compute_comb1(dt)
    comb2 <- compute_comb2(dt)
    output[[batch]] <- list(comb1, comb2)
  }
  saveRDS(output, sprintf("%s/combination.rds", path))
}





