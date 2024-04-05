set.seed(2024)
library(ggplot2)
library(tsibble)
library(dplyr)

for (data in c("mortality", "tourism")) {
  dt <- readRDS(sprintf("%s/data.rds", data))
  n <- NROW(dt$S)
  m <- NCOL(dt$S)
  dt_names <- colnames(dt$data)
  if (data == "tourism") {
    dt_select_middle <- sample(2:251, 1)
    dt_select_bottom <- (sample(1:76, 3) - 1) * 4 + 1:3
  } else {
    dt_select_middle <- sample(2:(n-m),1)
    dt_select_bottom <- sort(sample(1:m, 3))
  }
  mid_names <- colnames(dt$data[,dt_select_middle,drop=FALSE])
  bottom_names <- colnames(dt$data[,dt_select_bottom + n - m])
  dt_selected <- as.data.frame(dt$data[,c(1, dt_select_middle, dt_select_bottom  + n - m)])
  if (data == "tourism") {
    index <- make_yearmonth(year = rep(1998:2016, each = 12), month = rep(1:12, 19))
  } else {
    index <- make_yearmonth(year = rep(1999:2019, each = 12), month = rep(1:12, 21))
  }
  dt_selected <- dt_selected %>% 
    mutate(index = index) %>%
    tidyr::pivot_longer(cols = -index, names_to = "key", 
                        values_to = "y")
  
  p1 <- ggplot(data = dt_selected %>% filter(key == "Total") %>%
                 mutate(key = ifelse(data == "tourism", "Australian tourism", "U.S. Death"))) +
    geom_line(mapping = aes(x = index, y = y)) +
    ylab("") + xlab("") + 
    facet_wrap(~key)
  
  
  p2 <- ggplot(data = dt_selected %>% filter(key != "Total", key %in% mid_names)) +
    geom_line(mapping = aes(x = index, y = y)) +
    facet_wrap(~ key, ncol = 2, scales = "free") +
    ylab("") + xlab("") + guides(colour = "none")
  
  
  p3 <- gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
  pdf(sprintf("manuscript/figures/%s.pdf", data), width = 16 / 2.3 * 1.3, height=9 / 3 * 1.3)
  p4 <- ggplot(data = dt_selected %>% filter(key != "Total", key %in% bottom_names)) +
    geom_line(mapping = aes(x = index, y = y)) +
    facet_wrap(~ key, ncol = 3, scales = "free") +
    ylab("") + xlab("") + guides(colour = "none") +
    theme(panel.spacing.y = unit(1.5, "lines"))
  
  gridExtra::grid.arrange(p3, p4, nrow = 2, ncol=1, heights = c(1, 1))
  dev.off()
}



