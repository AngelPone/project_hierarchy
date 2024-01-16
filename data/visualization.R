set.seed(2024)
library(ggplot2)
library(tsibble)
library(dplyr)

for (data in c("mortality", "tourism")) {
  dt <- readRDS(sprintf("%s/data.rds", data))
  n <- NROW(dt$S)
  m <- NCOL(dt$S)
  dt_names <- readRDS(sprintf("%s/names.rds", data))
  dt_select_bottom <- sort(sample(1:m, 4))
  dt_selected <- as.data.frame(dt$data[,c(1, dt_select_bottom  + n - m)])
  colnames(dt_selected) <- c("Total", dt_names[dt_select_bottom])
  if (data == "tourism") {
    index <- make_yearmonth(year = rep(1998:2016, each = 12), month = rep(1:12, 19))
  } else {
    index <- make_yearmonth(year = rep(1999:2019, each = 12), month = rep(1:12, 21))
  }
  dt_selected <- dt_selected %>% 
    mutate(index = index) %>%
    tidyr::pivot_longer(cols = -index, names_to = "key", 
                        values_to = "y")
  
  pdf(sprintf("manuscript/figures/%s.pdf", data), width = 6, height=6)
  p1 <- ggplot(data = dt_selected %>% filter(key == "Total") %>%
                 mutate(key = ifelse(data == "tourism", "Australian tourism", "U.S. Death"))) +
    geom_line(mapping = aes(x = index, y = y)) +
    ylab("") + xlab("") + 
    facet_wrap(~key)
  p2 <- ggplot(data = dt_selected %>% filter(key != "Total")) +
    geom_line(mapping = aes(x = index, y = y)) +
    facet_wrap(~ key, ncol = 2, scales = "free") +
    ylab("") + xlab("") + guides(colour = "none")
  
  gridExtra::grid.arrange(p1, p2, nrow = 2, heights = c(1.5, 1.8))
  dev.off()
}



