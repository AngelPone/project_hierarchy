#' scripts to produce latex source code of table in manuscript
library(dplyr)

average_rmsse_table <- function(methods) {
  tourism <- read.csv("manuscript/figures/tourism/combination.csv", row.names = 1)
  mortality <- read.csv("manuscript/figures/mortality/combination.csv", row.names = 1)
  
  mortality <- mortality[mortality$method %in% methods,]
  tourism <- tourism[tourism$method %in% methods,]
  min_ <- c()
  for (i in 2:4) {
    min_ <- c(min_, tourism$method[which.min(tourism[,i])])
  }
  for (i in 2:4) {
    min_ <- c(min_, mortality$method[which.min(mortality[,i])])
  }
  rows <- list()
  top <- "\\begin{tabular}{lcccccc}\n\\toprule\n"
  head1 <- "Approach & \\multicolumn{3}{c}{tourism} & \\multicolumn{3}{c}{mortality} \\\\ \n"
  headrule <- "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n"
  head2 <- "& Top & Middle & Bottom & Top & Middle & Bottom \\\\ \\midrule\n"
  bottom <- "\n\\bottomrule\\end{tabular}"
  for (m in methods) {
    v_ <- c(as.numeric(tourism[tourism$method == m, 2:4]),
      as.numeric(mortality[mortality$method == m, 2:4]))
    bold_idx <- which(min_ == m)
    v_ <- sprintf("%.3f", v_)
    v_[bold_idx] <- paste0("\\textbf{", v_[bold_idx], "}")
    v_ <- sprintf("%s & %s \\\\", m, 
                  stringi::stri_join_list(as.list(v_), collapse=" & "))
    rows[[m]] <- v_
  }
  rows <- stringi::stri_join_list(rows, collapse = "\n")
  cat(top, head1, headrule, head2, rows, bottom)
}


methods1 <- c("Base", "Two-level", "Natural", 
              "TS-EUC-ME", "ER-EUC-ME", "TSF-EUC-ME", "ERF-EUC-ME",
              "TS-EUC-HC", "ER-EUC-HC", "ER-EUC-HC", "TSF-EUC-HC",
              "TS-DTW-ME", "TS-DTW-HC", "ER-DTW-ME", "ER-DTW-HC")

methods2 <- c("Base", "Two-level", "Natural", "TS-DTW-HC", "TSF-EUC-HC", "Combination", "Stack")

average_rmsse_table(methods1)
average_rmsse_table(methods2)
