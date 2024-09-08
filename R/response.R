# zero ratio
batch_lengths <- list(tourism=120, mortality=144)
for (path in c("tourism", "mortality")) {
  batch_length <- batch_lengths[[path]]
  n0 <- 0
  n <- 0
  for (i in 0:batch_length) {
    dt <- readRDS(sprintf("%s/batch_%d.rds", path, i))
    natural <- Filter(\(x) x$cluster=="natural", dt$nl)[[1]]
    n0 <- n0 + sum(natural$rf < 0)
    n <- n + length(natural$rf)
  }
  print(sprintf("%s percent of 0 in %s dataset", n0 / n, path))
}
