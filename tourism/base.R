# configuration
args <- commandArgs(trailingOnly = TRUE)
batch <- as.integer(args[[1]])
path <- "tourism"
bfmethod <- "arima"
num.cores <- 8
set.seed(42)

cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

source("R/construct_hierarchy.R", chdir = T)


dt <- readRDS("tourism/data.rds")

## load dataset
data <- hts(rbind(rep(1, 304), diag(304)), 
            bts = dt$data[1:(216-batch), 252:555],
            tts = dt$data[(217-batch):(228-batch), 252:555])

# Workflow

## load base forecast store
BASEFORECAST_STORE <- list()
BASEFORECAST_STORE[[bfmethod]] <- list()
output <- create_new_output()


## bottom and total forecast
data <- forecast(data, bfmethod)
accs <- evaluate.hts(data, metrics, type = "base")
output <- add_result(output, "", "", "", accs)
data <- reconcile.all(data)
accs <- evaluate.hts(data, metrics, type = "rf")
output <- add_result(output, "", "", "", accs)
FEATURES <- NULL
saveResult()


