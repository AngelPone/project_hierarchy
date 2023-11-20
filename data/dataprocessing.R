library(dplyr)
rm(list = ls())
dt <- read.csv("data/mortality_1999-2019.csv") %>%
  select(-Notes, -X, -Year.Code, -Month, -Population, -Crude.Rate) %>%
  filter(ICD.10.113.Cause.List.Code != "") %>%
  rename(year=Year, month=Month.Code, desc=ICD.10.113.Cause.List,
         code=ICD.10.113.Cause.List.Code, deaths=Deaths)


dt_total <- read.delim("data/mortality_1999-2019_total.txt") %>%
  select(-Notes, -Year.Code, -Month, -Population, -Crude.Rate) %>%
  rename(year=Year, month=Month.Code, deaths=Deaths) %>%
  filter(year <= 2019)

ICDcodes <- unique(dt$desc)
ICDcodes <- ICDcodes[ ! (ICDcodes %in% c("#COVID-19 (U07.1)", "#Enterocolitis due to Clostridium difficile (A04.7)"))]
ICDcodes <- data.frame(ICDcodes, level = c(
  rep(1, 4), 2, 2, rep(1, 13), rep(2, 18), rep(3, 5), 2, rep(1, 4), 2, 2, rep(1, 4),
  2, rep(3, 4), rep(4, 3), 5, 5, 3, rep(4, 4), rep(2, 4), 3, 3, 1, 1, 2, 2, 1, 2, 2,
  1, rep(2, 4), rep(1, 7), 2, 2, 1, 1, rep(2, 4), rep(1, 4), 2, 2, rep(1, 5), 2, rep(3,3),
  2, rep(3, 6), 1, 2, 2, 1, 2, 2, 1, 1, 2, 2, rep(1, 2)
))


ICD_leaf <- vector(length = NROW(ICDcodes))
for (i in 1:(NROW(ICDcodes) - 1)) {
  if (as.integer(ICDcodes[i, 2]) >= as.integer(ICDcodes[i+1, 2])) {
    ICD_leaf[i] <- TRUE
  }
}
ICD_leaf[NROW(ICDcodes)] <- TRUE

ICDcodes$leaf <- ICD_leaf
rm(ICD_leaf)

m <- sum(ICDcodes$leaf)
n <- NROW(ICDcodes)
# summing matrix
S <- matrix(0, NROW(ICDcodes), m)
bottom_idx <- NULL
idx <- 0
for (i in 1:NROW(S)) {
  if (ICDcodes$leaf[i]) {
    idx <- idx + 1
    bottom_idx <- c(bottom_idx, idx)
  } else {
    bottom_idx <- c(bottom_idx, 0)
  }
}
ICDcodes$bottom_idx <- bottom_idx

find_leaves <- function(start_idx) {
  leave_idx <- c()
  level <- ICDcodes$level[start_idx]
  while(TRUE) {
    start_idx <- start_idx + 1
    if (ICDcodes$level[start_idx] <= level) {
      break
    } else if (ICDcodes$leaf[start_idx]){
      leave_idx <- c(leave_idx, ICDcodes$bottom_idx[start_idx])
    }
  }
  return(leave_idx)
}


for (i in 1:NROW(S)) {
  if (ICDcodes$leaf[i]) {
    S[i, ICDcodes$bottom_idx[i]] <- 1
  } else {
    S[i, find_leaves(i)] <- 1
  }
}



print_leaf <- function(i) {
  print("parent Node:")
  print(ICDcodes$ICDcodes[which(!ICDcodes$leaf)[i]])
  print("Nodes:")
  bottom_idx <- which(S[which(!ICDcodes$leaf)[i],] == 1)
  
  print(ICDcodes$ICDcodes[which(ICDcodes$bottom_idx %in% bottom_idx)])
}
# check hierarchy
print_leaf(2)

rownames(S) <- ICDcodes$ICDcodes
colnames(S) <- ICDcodes$ICDcodes[ICDcodes$bottom_idx>0]


combine_cols <- dt %>% group_by(desc) %>% summarise(n_suppressed = sum(deaths == "Suppressed")) %>%
  filter(n_suppressed > 0) %>% left_join(ICDcodes %>% rename(desc = ICDcodes))

find_leaves <- function(name) {
  S_row <- S[name,]
  if (sum(S_row) == 1) {
    return(NULL)
  } else {
    return(colnames(S)[which(S_row == 1)])
  }
}

allSeries <- combine_cols$desc
series <- allSeries[combine_cols$level == 1]
# remove from S

rowsToRemove <- NULL
colsToRemove <- NULL
for (s in series) {
  leves <- find_leaves(s)
  rowsToRemove <- c(rowsToRemove, s, leves)
  if (length(leves) > 0) {
    colsToRemove <- c(colsToRemove, leves)
  } else {
    colsToRemove <- c(colsToRemove, s)
  }
}
rowsToRemoveIdx <- which(rownames(S) %in% rowsToRemove)
colsToRemoveIdx <- which(colnames(S) %in% colsToRemove)
S <- S[-rowsToRemoveIdx, -colsToRemoveIdx]
S <- rbind(S, vector("numeric", NCOL(S)))
rownames(S)[NROW(S)] <- "Other"
S <- cbind(S, vector("numeric", NROW(S)))
colnames(S)[NCOL(S)] <- "Other"
S[NROW(S), NCOL(S)] <- 1


# remove from dt
otherSeries <- setdiff(ICDcodes$ICDcodes[ICDcodes$level == 1], rowsToRemove)
deathsOther <- dt %>% filter(desc %in% otherSeries) %>%
  group_by(year, month) %>% summarise(deaths = sum(as.integer(deaths))) %>%
  left_join(dt_total %>% rename(deaths_total=deaths), by = c("year", "month")) %>%
  mutate(death_Other = deaths_total - deaths)
dt <- dt %>% filter(!(desc %in% rowsToRemove))
dt <- dt %>% rbind(deathsOther %>% select(year, month, deaths=death_Other) %>%
                     mutate(desc = "Other", code = "Other"))

ICDcodes <- ICDcodes %>% filter(!(ICDcodes %in% rowsToRemove))
combine_cols <- combine_cols %>% filter(!(desc %in% rowsToRemove))


find_parent <- function(name) {
  alts <- rownames(S)[which(S[, name] == 1)]
  for (alt in alts) {
    if (ICDcodes$level[ICDcodes$ICDcodes == alt] == (ICDcodes$level[ICDcodes$ICDcodes == name]-1)) {
      return(alt)
    }
  }
}

# others

allSeries <- combine_cols$desc
for (code in combine_cols$desc) {
  if (!(code %in% allSeries)) { next }
  parent <- find_parent(code) 
  print(sprintf("current series is %s", code))
  print(paste0("parent is ", parent))
  allLeaves <- find_leaves(parent)
  otherLeaves <- setdiff(allLeaves, allSeries)
  combineLeaves <- intersect(allLeaves, allSeries)
  print(paste0("series to remove is ", combineLeaves))
  print(paste0("series to minus is ", otherLeaves))
  # 
  allSeries <- allSeries[!(allSeries %in% combineLeaves)]
  ICDcodes <- ICDcodes %>% filter(!(ICDcodes %in% combineLeaves))
  combine_cols <- combine_cols %>% filter(!(desc %in% combineLeaves))
  
  if (length(combineLeaves) == 1) {
    dt$deaths[dt$desc == combineLeaves] <-
      dt %>% filter(desc == parent) %>% mutate(deaths = as.integer(deaths)) %>%
      left_join(dt %>% filter(desc %in% otherLeaves) %>% group_by(year, month) %>%
                  summarise(deaths = sum(as.integer(deaths))), by = c("year", "month")) %>%
      mutate(deaths = deaths.x-deaths.y) %>% pull(deaths)
    next
  }
  
  
  # S
  newname <- paste0("Combine-", parent)
  rowsToRemoveIdx <- which(rownames(S) %in% combineLeaves)
  colsToRemoveIdx <- which(colnames(S) %in% combineLeaves)
  S <- S[-rowsToRemoveIdx, -colsToRemoveIdx]
  S <- rbind(S, vector("numeric", NCOL(S)))
  rownames(S)[NROW(S)] <- newname
  S <- cbind(S, vector("numeric", NROW(S)))
  colnames(S)[NCOL(S)] <- newname
  S[parent, newname] <- 1
  S[newname, newname] <- 1
  
  # dt
  deathsOther <- dt %>% filter(desc %in% otherLeaves) %>%
    mutate(deaths = as.integer(deaths)) %>%
    group_by(year, month) %>% summarise(deaths = sum(deaths)) %>%
    left_join(dt %>% filter(desc == parent) %>% 
                rename(death_total = deaths) %>% 
                mutate(death_total = as.integer(death_total)), by = c("year", "month")) %>%
    mutate(death_other = death_total - deaths, code = newname,
           desc = newname) %>%
    select(year, month, desc, code, deaths=death_other)
  
  dt <- dt %>% filter(!(desc %in% combineLeaves)) %>%
    rbind(deathsOther)
}


dt <- dt %>% mutate(deaths = as.integer(deaths))

dt <- dt %>% filter(!(desc %in% c("#Enterocolitis due to Clostridium difficile (A04.7)", 
                                  "#COVID-19 (U07.1)")) )

# test summation equality
for (code in unique(dt$desc)) {
  
  allLeaves <- find_leaves(code)
  
  if (length(allLeaves) > 0) {
    diff <- dt %>% filter(desc %in% allLeaves) %>%
      group_by(year, month) %>% summarise(deaths_sum = sum(deaths)) %>%
      left_join(dt %>% filter(desc == .env$code), by = c("year", "month")) %>%
      mutate(diff = deaths - deaths_sum) %>%
      pull(diff)
    stopifnot( all( diff == 0))
  }
}


# test summing matrix

wide_dt <- dt %>% select(month, desc, deaths) %>%
  tidyr::pivot_wider(id_cols = month, names_from = desc, values_from = deaths)

dot <- as.matrix(wide_dt[,colnames(S)]) %*% t(S)

all(as.matrix(wide_dt[colnames(dot)]) == dot)


output <- wide_dt[colnames(S)]
colnames(output) <- data.frame(desc = colnames(S)) %>% left_join(dt %>% select(desc, code) %>% unique(), by = "desc") %>%
  pull(code)


# remove all series whose all values are zero
# GR113-012
output <- output[,-which(colMeans(output) == 0)]
colnames(output)[NCOL(output)] <- "Combine-GR113-097"

S <- S[-which(rownames(S) == "#Acute poliomyelitis (A80)"),-which(colnames(S) == "#Acute poliomyelitis (A80)")]

colnames(S) <- data.frame(desc = colnames(S)) %>% 
  left_join(dt %>% select(desc, code) %>% unique(), by = "desc") %>%
  pull(code)

rownames(S) <- data.frame(desc = rownames(S)) %>% 
  left_join(dt %>% select(desc, code) %>% unique(), by = "desc") %>%
  pull(code)

colnames(S)[NCOL(S)] <- "Combine-GR113-097"
rownames(S)[NROW(S)] <- "Combine-GR113-097"


S <- rbind(S[!(rownames(S) %in% colnames(S)),],
           S[colnames(S),])
S <- rbind(rep(1, 98), S)
rownames(S)[1] <- "Total"
output <- as.matrix(output) %*% (t(S))
saveRDS(list(S=S, data=output), "mortality/data.rds")



