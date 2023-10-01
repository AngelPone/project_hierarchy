batch <- 0
library(dplyr)
library(ggplot2)
mortality <- readRDS(sprintf("mortality/store_%s.rds", batch))$data
tourism <- readRDS(sprintf("tourism/store_%s.rds", batch))$data

visual <- function(dt, idx) {
  
  resid <- unname(as.matrix(dt$resid[,idx+1]))
  ts <- unname(as.matrix(dt$bts[,idx]))
  
  resid <- apply(resid, 2, function(x){
    (x - mean(x))/sd(x)
  })
  ts <- apply(ts, 2, function(x){
    (x-mean(x))/sd(x)
  })
  
  colnames(ts) <- paste("Series", 1:NCOL(ts))
  colnames(resid) <- paste("Series", 1:NCOL(ts))
  
  data.frame(resid) %>% mutate(type = "error", x=1:NROW(resid)) %>%
    rbind(data.frame(ts) %>% mutate(type = "ts", x=1:NROW(resid))) %>%
    tidyr::pivot_longer(cols = -c("type", "x")) %>%
    ggplot(mapping = aes(x=x, y=value)) +
    geom_line() +
    facet_wrap(~ type + name, nrow = 2)
}
visual(tourism, sample(304, 4))
visual(mortality, sample(98, 4))

