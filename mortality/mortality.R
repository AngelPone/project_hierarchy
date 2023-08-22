# configuration
args <- commandArgs(trailingOnly = TRUE)
batch <- as.integer(args[[1]])
path <- "mortality"
bfmethod <- "arima"
num.cores <- 8
cl <- parallel::makeCluster(num.cores)
doParallel::registerDoParallel(cl)

source("R/construct_hierarchy.R", chdir = T)


mortality <- read.csv("data/mortality.csv")
colnames(mortality) <- stringi::stri_replace_all_fixed(colnames(mortality), ".", "-")


## load dataset
data <- hts(rbind(rep(1, 98), diag(98)), 
            bts = unname(as.matrix(mortality))[1:(240-batch),],
            tts = unname(as.matrix(mortality))[(241 - batch):(252 - batch),])

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

saveResult()


