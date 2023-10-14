batch <- 0

data <- readRDS(sprintf("tourism/ets/store_%s.rds", batch))$distance

objs <- list()

REPRESENTORS <- c("ts", "error", "ts.features", "error.features", "forecast")
DISTANCES <- c("euclidean", "dtw", "negcor", "cor", "uncorrelation")

for (representor in REPRESENTORS) {
  objs[[representor]] <- list()
  for (distance in DISTANCES) {
    objs[[representor]][[distance]] <- foreach(n_clusters = 1:15, .packages = "cluster") %do% {
      
      distance_mat <- data[[representor]][[distance]]
      pam(distance_mat, k = n_clusters, diss=TRUE)$objective
    }
  }
}


for (representor in REPRESENTORS) {
  for (distance in DISTANCES) {
    objs[[representor]][[distance]] <- do.call(rbind, objs[[representor]][[distance]])
  }
}


plot_dt <- NULL
for (representor in REPRESENTORS) {
  x <- objs[[representor]][[DISTANCES[1]]]
  x <- tibble(values=t(apply(x, 1, function(g){g/x[1,]}))[,2], ncluster=1:15, label=DISTANCES[1])
  for (distance in DISTANCES[2:length(DISTANCES)]) {
    x2 <- objs[[representor]][[distance]]
    x2 <- tibble(values = t(apply(x2, 1, function(g){g/x2[1,]}))[,2],
                 ncluster = 1:15, label = distance)
    x <- rbind(x, x2)
  }
  plot_dt <- rbind(plot_dt, x %>% mutate(representor = .env$representor))
}


library(ggplot2)
ggplot(plot_dt) +
  geom_line(aes(x=ncluster, y=values, color = label, group = label)) +
  facet_wrap(~representor)



