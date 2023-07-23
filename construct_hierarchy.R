

#' hts
#' 
#' S: summing matrix
#' bts: bottom level series
#' nl: new level
#" basef: base forecasts
#' rf: reconciled forecasts
hts <- function(S, bts) {
  structure(
    list(
      S = S, bts = bts,
      nl = list(),
      basef = NULL,
      rf = NULL
    ),
    class = "hts"
  )
}



#' function to forecast hts
forecast.hts <- function(x, f){
  if (is.null(x$basef)) {
    all_ts <- x$bts %*% t(x$S)
    x$basef <- apply(all_ts, 2, f)
  }
  if (!is.null(hts$nl$S) && is.null(hts$nl$basef)) {
    all_nts <-x$bts %*% t(x$nl$S)
    x$nl$basef <- apply(all_ts, 2, f)
  }
}

#' add new level
updateLevel.hts <- function(x, new_S, meta) {
  append(x$nl, list(S = S, basef=NULL, rf=NULL, rf_nl=NULL, meta=meta))
}

evaluate.hts <- function(x, metric = "RMSE") {
  if (is.null(x$rf)) {
    stop("Reconciled forecasts are null")
  }
  
  if (metric == "RMSE") {
    for (level in x$nl) {
      if (is.null(level$rf)) next
      level$acc <- sqrt(rowMeans((level$rf - x$bts)^2))
    }
  }
}


#' Euclidean distance
distance.euclidean <- function(x, y) { (x-y)^2 }



#' construct hierarchy
#' 
#' @param cluster method of clustering
build_level <- function(
    hts,
    distance,
    by = c("error", "tsfeat", "ts", "errfeat"),
    cluster = c("km", "hclut")) {
  
  by <- match.arg(by)
  cluster <- match.arg(cluster)
  
  if (by == "ts") {
    
  }
  
}