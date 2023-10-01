library(dplyr)
library(tsibble)

prison_orgi <- readr::read_csv("https://OTexts.com/fpp3/extrafiles/prison_population.csv") |>
  mutate(Quarter = yearquarter(Date)) |>
  select(-Date)  |>
  as_tsibble(key = c(Gender, Legal, State, Indigenous),
             index = Quarter) |>
  relocate(Quarter)


prison <- prison_orgi %>% group_by(State, Gender) %>% summarise_at("Count", sum)

States <- unique(prison$State)
Genders <- unique(prison$Gender)

prison <- prison %>% tibble() %>%
  mutate(id_col = paste0(State, "-", Gender)) %>%
  select(-c("State", "Gender")) %>%
  tidyr::pivot_wider(names_from = id_col, values_from = "Count", id_cols = Quarter)

m <- length(States) * length(Genders)
S <- matrix(0, 1 + m + length(States) + length(Genders), m)
S[1,] <- rep(1, length(m))

bottom_names <- colnames(prison)[2:NCOL(prison)]
bottom_states <- sapply(strsplit(bottom_names, "-"), function(x){ x[[1]] })
bottom_genders <- sapply(strsplit(bottom_names, "-"), function(x){ x[[2]] })

for ( i in 1:length(States)) {
  S[i+1, which(bottom_states == States[i])] <- 1
}

for (i in 1:length(Genders)) {
  S[(1+length(States)+i), which(bottom_genders == Genders[i])] <- 1
}

S[(1+length(States) + length(Genders) + 1) : NROW(S),] <- diag(m)


series_names <- c("Total", States, Genders, bottom_names)
bts <- unname(as.matrix(prison[2:NCOL(prison)]))


allts <- bts %*% t(S)


# assert summing matrix is correct
all(sapply(States, function(state){
  all.equal(prison_orgi %>% filter(State == state) %>% group_by() %>% summarise_at("Count", sum) %>%
              pull(Count), allts[, which(series_names == state)])
}))


all(sapply(Genders, function(gender){
  all.equal(prison_orgi %>% filter(Gender == gender) %>% group_by() %>% summarise_at("Count", sum) %>%
              pull(Count), allts[, which(series_names == gender)])
}))

all.equal(prison_orgi %>% group_by() %>% summarise_at("Count", sum) %>% pull(Count), allts[,1])

# save Data
source("R/construct_hierarchy.R", chdir = T)

num.cores <- 8

cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

metrics <- c("rmse", "mae")


write.csv(S, "prison/S.csv", row.names = FALSE)

for (i in 0:7) {
  
  BASEFORECAST_STORE <- list()
  BASEFORECAST_STORE[["arima"]] <- list()
  FEATURES <- NULL
  
  store_path <- sprintf("prison/store_%s.rds", i)
  
  train <- bts[1:(NROW(bts)-i-4),,drop=FALSE]
  test <- bts[(NROW(bts)-i-3):(NROW(bts) - i), ,drop=FALSE]
  data <- hts(S[c(1, (NROW(S) - NCOL(S) + 1):NROW(S)), ], train, test)
  data <- forecast(data, "arima", h=NROW(test), frequency=4)
  
  output <- create_new_output()
  
  accs <- evaluate.hts(data, metrics, type = "base")
  output <- add_result(output, "", "", "", accs)
  data <- reconcile.all(data)
  accs <- evaluate.hts(data, metrics, type = "rf")
  output <- add_result(output, "", "", "", accs)
  
  saveResult()
}







