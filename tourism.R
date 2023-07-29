source("construct_hierarchy.R")
library(forecast)

dt <- read.csv("~/Documents/projects/datasets/TourismLarge/data.csv")
series_names <- colnames(dt)
dt <- dt[serise_names[endsWith(serise_names, "All") & (stringi::stri_length(serise_names) == 6)]]
rm(series_names)

dt <- hts(S = rbind(rep(1, NCOL(dt)), diag(NCOL(dt))),
          bts = unname(as.matrix(dt)))


# ts, Euclidean distance, kmeans


dt_1 <- build_level(dt, 
            representator = representator.ts, 
            distance = distance.euclidean,
            cluster = cluster.kmedoids,
            n_clusters = c(10, 15))

dt_1 <- forecast(dt_1, function(x){
  forecast(auto.arima(x), h=12)$mean
})
