# Check that other variables are constant
check_unused_variables <- function(data, usednames) {
  for(name in setdiff(names(pairdata), usednames)) {
    if(!all(pairdata[[name]] == pairdata[[name]][1]))
      stop(paste("Don't know what to do with '", name, "' variable.."))
  }
}

check_activity_split <- function(...) {
  acts <- list(...)
  vars <- c(acts[[1]]$factors, acts[[1]]$by)

  factorcombinations <- do.call(rbind, lapply(acts, function(x) {
    d <- x$transmat
    d$varlevels <- NULL
    d$transmat <- NULL
    d
  }))
  dups <- duplicated(factorcombinations)
  dups2 <- duplicated(factorcombinations[vars])
  stopifnot(identical(dups, dups2))
  if(any(dups)) {
    print("Duplicate activies found for:")
    print(factorcombinations[dups,])
  }
}
