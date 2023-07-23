library()

evaluate <- function(obj){
  rf <- obj$rf
  S <- rbind(rep(1, 4), diag(4))
  y_agg <- (obj$y %*% t(S))[(TRAIN_WINDOW+1):TIME_WINDOW,]
  tmp <- sapply(rf, function(y){
    tmp <- sapply(y, function(k){
      sqrt(colMeans((k - y_agg)^2))
    }, simplify = "array")
    colnames(tmp) <- names(y)
    rownames(tmp) <- paste0('s', 1:5)
    tmp
  }, simplify = "array")
  dimnames(tmp)[[3]] <- c("S1", "S2", "S3")
  tmp
}

foo <- evaluate(objs[[1]])

evaluate_basef <- function(obj){
  S <- rbind(rep(1, 4), diag(4))
  y_agg <- (obj$y %*% t(S))[(TRAIN_WINDOW+1):TIME_WINDOW,]
  f <- cbind(obj$f$total, obj$f$bottom)[(TRAIN_WINDOW+1):TIME_WINDOW,]
  colMeans((y_agg - f)^2)
}


objs <- readRDS("data/objs.rds")

accs <- sapply(objs, evaluate, simplify = "array")
accs_basef <- sapply(objs, evaluate_basef, simplify = "array")

# saveRDS(accs, 'data/accuracy.rds')

library(tsutils)

series <- paste0("s", 1:5)
methods <- c("ols", "wlss", "wlsv", "shrinkage")
groups <- c("S1", "S2", "S3")

pdf("mean.pdf", height = 10, width = 10)
par(mfrow = c(2, 2))
for (method in methods){
  f1 <- cbind(t(apply(accs, c(2,3,4), mean)[method,,]), apply(accs_basef, 2, mean))
  colnames(f1)[4] <- "base"
  nemenyi(f1, plottype = "vmcb", main = sprintf("MCB Test on mean RMSE of the hierarchy for method %s", method))
}
dev.off()

series_mcb <- function(acc, S = 1){
  sname <- paste0("S", S)
  pdf(paste0(sname, ".pdf"), height = 10, width = 10)
  par(mfrow = c(2, 2))
  for (method in methods){
    f1 <- cbind(t(accs[S,method,,]), accs_basef[S,])
    colnames(f1)[4] <- "base"
    nemenyi(f1, plottype = "vmcb", main = sprintf("MCB Test on mean RMSE of Series %d for method %s", S, method))
  }
  dev.off()
}


series_mcb(accs, 1)
series_mcb(accs, 2)
series_mcb(accs, 3)
series_mcb(accs, 4)
series_mcb(accs, 5)


accs_summary <- apply(accs, c(1,2,3), mean)

output_tbs <- list()
for (method in methods){
  output_tbs[[method]] <- rbind(accs_summary[,method,], colMeans(accs_summary[,method,]))
  rownames(output_tbs[[method]])[6] <- "average"
  print(method)
  print(output_tbs[[method]])
}






