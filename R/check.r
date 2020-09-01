# Check that other variables are constant
check_unused_variables <- function(data, usednames) {
  for(name in setdiff(names(data), usednames)) {
    if(!all(data[[name]] == data[[name]][1]))
      stop(paste("Don't know what to do with '", name, "' variable.."))
  }
}

check_activity_split <- function(acts, ...) {
  if(length(list(...))) acts <- c(list(acts), list(...))
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
    print(paste0("Duplicate entries found for activity '", acts[[1]]$actprobname, "':"))
    print(factorcombinations[dups,])
    return(TRUE)
  }
  FALSE
}
