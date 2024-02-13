rm(list=ls())
dt <- readRDS("simulation/simulation.rds")
source("R/metrics.R")

rmsse <- function(pred, obs, hist) {
  obs <- cbind(rowSums(obs), obs)
  hist <- cbind(rowSums(hist), hist)
  mean(sapply(1:NCOL(pred), function(x){
    metric.rmsse(pred[,x], obs[,x], hist[,x])
  }))
}


no_hierarchy <- 1
best_hierarchy <- 2
permute_6 <- 3:102
permute_3 <- 103:202
permute_2 <- 203:302
trend_dir <- 303
season <- 304
trned_exis <- 305


compare_random <- function(idx_orig, idx_random) {
  f <- dt$acc
  all_rmsse <- NULL
  for (i in seq_along(dt$acc)) {
    tts <- t(dt$series[[i]][,133:144])
    bts <- t(dt$series[[i]][,1:132])
    f_orig_ <- f[[i]][[idx_orig]]
    f_permu_ <- f[[i]][idx_random]
    
    rmsse_orig <- rmsse(f_orig_, tts, bts)
    rmsse_permu_ <- sapply(iterators::iter(f_permu_), function(g){
      rmsse(g, tts, bts)
    })
    all_rmsse <- rbind(all_rmsse, c(rmsse_orig, rmsse_permu_))
  }
  all_rmsse
}

test <- function(mat, name) {
  for (i in 1:NROW(mat)) {
    mat[i,2:101] <- sort(mat[i,2:101])
  }
  colnames(mat) <- c(name, 1:100)
  test_ <- tsutils::nemenyi(mat, plottype = "none")
  which(colnames(test_$interval) == name)
}

# natural hierarchy vs its counterpart
natural_ <- compare_random(best_hierarchy, permute_6)
test(natural_, "natural")
# 39

# season
season_ <- compare_random(season, permute_2)
test(season_, "season")
# 53

# trend existence
trned_exis_ <- compare_random(trned_exis, permute_2)
test(trned_exis_, "trend existence")
# 48

# trend direction
trend_dir_ <- compare_random(trend_dir, permute_3)
test(trend_dir_, "trend direction")
# 27


method_idx <- c(2, 303, 304, 305)
all_rmsse <- NULL
for (i in seq_along(dt$acc)) {
  print(i)
  tts <- t(dt$series[[i]][,133:144])
  bts <- t(dt$series[[i]][,1:132])
  f_orig_ <- apply(simplify2array(dt$acc[[i]][method_idx]), c(1,2), mean)
  f_permu_ <- lapply(1:100, function(j) {
    apply(simplify2array(
      dt$acc[[i]][c(permute_6[j], permute_3[j], permute_2[j], permute_2[j])]
    ), c(1,2), mean)
  })
  
  rmsse_orig <- rmsse(f_orig_, tts, bts)
  rmsse_permu_ <- sapply(iterators::iter(f_permu_), function(g){
    rmsse(g, tts, bts)
  })
  all_rmsse <- rbind(all_rmsse, c(rmsse_orig, rmsse_permu_))
}

# 51
test(all_rmsse, "average")

for (i in 1:500) {
  all_rmsse[,2:101] <- sort(all_rmsse[,2:101])
}

mat <- cbind(natural_[,1], season_[,1], trend_dir_[,1], trned_exis_[,1], all_rmsse)
colnames(mat) <- c("Natural", "Season", "Trend direction", "trend existence", "average", 1:100)
test_ <- nemenyi(mat, plottype = "none")

which(colnames(test_$intervals) == "Natural") #51
which(colnames(test_$intervals) == "Season") # 56
which(colnames(test_$intervals) == "Trend direction") #49
which(colnames(test_$intervals) == "trend existence") #55
which(colnames(test_$intervals) == "average") #52
