library(FoReco)
library(forecast)
library(dplyr)
source("R/reconciliation.R")

mortality <- readRDS("tourism/store_5.rds")
bts <- mortality$data$bts
tts <- mortality$data$tts
S <- mortality$data$S

# check the base forecasts

check_bf <- function() {
  idx <- sample(length(mortality$bfstore$ets), 1)
  series_name <- names(mortality$bfstore$ets)[idx]
  S <- as.numeric(strsplit(series_name, "")[[1]])
  series <- bts %*% matrix(S, ncol=1)
  series <- ts(series[,1], frequency = 12)
  tts <- ts((tts %*% matrix(S, ncol=1))[,1], frequency = 12, start = c(20, 8))
  ets_mdl <- ets(series)
  
  # check forecasts
  fcasts <- forecast(ets_mdl, h=12)
  print(sprintf("base forecasts: %s", all.equal(fcasts$mean, mortality$bfstore$ets[[series_name]]$basef)))
  resid <- residuals(ets_mdl)
  print(sprintf("residuals: %s", all.equal(resid, mortality$bfstore$ets[[series_name]]$resid)))
  
  # check plots
  plot(series, main = series_name)
  lines(fcasts$mean, col="red")
  lines(fitted(ets_mdl), col="blue")
  lines(tts, col="green")
}

check_bf()

# check reconciliation
check_rf <- function(idx = NULL) {
  if (is.null(idx)) {
    idx <- sample(length(mortality$output$cluster), 1)
  }
  m <- NCOL(mortality$data$bts)
  S <- mortality$output$other[[idx]]$S
  if (!is.null(S)) {
    S <- as.matrix(S)
  }
  C <- rbind(rep(1, m), S)
  S <- rbind(C, diag(m))
  tts <- mortality$data$tts %*% t(S)
  
  base <- lapply(iterators::iter(S, by="row"), function(x){
    series_name <- do.call(paste0, as.list(x))
    mortality$bfstore$ets[[series_name]]$basef
  }) %>% do.call(cbind, .)
  
  res <- lapply(iterators::iter(S, by="row"), function(x){
    series_name <- do.call(paste0, as.list(x))
    mortality$bfstore$ets[[series_name]]$resid
  }) %>% do.call(cbind, .)
  
  ols_foRec <- htsrec(base, C=C, comb = "ols")$recf
  ols_my <- reconcile.ols(S, base, res)
  print(sprintf("OLS difference: %s", max(abs(ols_my - ols_foRec))))

  wlss_foRec <- htsrec(base, C=C, comb = "struc")$recf
  wlss_my <- reconcile.wlss(S, base, res)
  print(sprintf("Structural Scaling difference: %s", max(abs(wlss_my - wlss_foRec))))

  wlsv_foRec <- htsrec(base, C=C, comb = "wls", res = res)$recf
  wlsv_my <- reconcile.wlsv(S, base, res)
  print(sprintf("Variance Scaling difference: %s", max(abs(wlsv_my - wlsv_foRec))))
  
  mint_foRec <- htsrec(base, C=C, comb = "shr", res = res)$recf
  mint_my <- reconcile.mint(S, base, res)
  print(sprintf("Shrinkage difference: %s", max(abs(mint_my - mint_foRec))))
  
  # accuracy check
  series_idx <- c(1, (NROW(C)+1):NROW(S))
  
  for (rf_method in c("ols", "wlsv", "wlss", "mint")) {
    orig_acc <- mortality$output$accuracy[[idx]]$rmse[[rf_method]]
    orig_acc <- c(orig_acc$total, orig_acc$bottom)
    new_acc <- sapply(series_idx, function(x) {
      sqrt(mean((get(paste0(rf_method, "_foRec"))[, x] - tts[, x])^2))
    })
    print(sprintf("%s accuracy difference: %s", rf_method, max(abs(orig_acc - new_acc))))
  }
}

check_rf()


check_relative <- function() {
  no_hierarchy_idx <- 2
  natural_idx <- which(mortality$output$cluster == "natural")
  
  for (rf_method in c("ols", "wlss", "wlsv", "mint")) {
    no_hierarchy <- mortality$output$accuracy[[no_hierarchy_idx]]$rmse[[rf_method]]
    natural <- mortality$output$accuracy[[natural_idx]]$rmse[[rf_method]]
    
    no_hierarchy <- c(no_hierarchy$total, no_hierarchy$bottom)
    natural <- c(natural$total, natural$bottom)
    
    print(sprintf("%s, natural better than no hierarchy: %s/%s", rf_method, sum(natural < no_hierarchy), length(natural)))
  }
  
  for (rf_method in c("ols", "wlss", "wlsv", "mint")) {
    no_hierarchy <- mortality$output$accuracy[[no_hierarchy_idx]]$rmse[[rf_method]]
    natural <- mortality$output$accuracy[[natural_idx]]$rmse[[rf_method]]
    
    no_hierarchy <- c(no_hierarchy$total, mean(no_hierarchy$bottom))
    natural <- c(natural$total, mean(natural$bottom))
    
    print(sprintf("%s, total average performance: natural: %s no_hierarchy: %s", rf_method, natural[1], no_hierarchy[1]))
    print(sprintf("%s, bottom average performance: natural: %s no_hierarchy: %s", rf_method, natural[2], no_hierarchy[2]))
  }
}
check_relative()



