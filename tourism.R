source("construct_hierarchy.R")
source("representator.R")
source("distance.R")
source("clustering.R")
source("reconciliation.R")

library(forecast)

dt <- read.csv("~/Documents/projects/datasets/TourismLarge/data.csv")
series_names <- colnames(dt)
dt <- dt[series_names[endsWith(series_names, "All") & (stringi::stri_length(series_names) == 6)]]
rm(series_names)

dt <- hts(S = rbind(rep(1, NCOL(dt)), diag(NCOL(dt))),
          bts = unname(as.matrix(dt))[1:216,],
          tts = unname(as.matrix(dt))[217:228,])
dt <- forecast(dt, f.arima)
dt <- reconcile.wlss(dt)

acc_base <- evaluate.hts(dt)

# ts, Euclidean distance, kmeans
dt_1 <- build_level(dt, 
            representator = representator.ts, 
            distance = distance.euclidean,
            cluster = cluster.kmedoids,
            n_clusters = c(10))

dt_1 <- forecast(dt_1, f.arima)
dt_1 <- reconcile.wlss(dt_1)

acc_g1 <- evaluate.hts(dt_1)



dt_2 <- build_level(dt_1, 
                    representator = representator.ts, 
                    distance = distance.euclidean,
                    cluster = cluster.kmedoids,
                    n_clusters = c(15), keep_old = TRUE)
dt_2 <- forecast(dt_2, f.arima)
dt_2 <- reconcile.wlss(dt_2)
acc_g2 <- evaluate.hts(dt_2)


View(cbind(acc_base, acc_g1, acc_g2))
