#' Define activity
#'
#' @param name Name of activity used in reporting
#' @param statespace Names of variables where changes happen
#' @param probname Name of activity in activity probabilities data
#' @examples
#' define_activity("nomanagement", c("vol", "age"))
#' @export
define_activity <- function(name, statespace, probname=name) {
  list(actname=name, actprobname=probname, statespace0=paste0(statespace, "0"), statespace1=paste0(statespace, "1"))
}

#' Add a statespace to an activity
#'
#'
#' Possible states are found using statespace0 and statespace1.
#' Possible factor combinations are found from data. Thus the data should contain at least one observation for each
#' interesting combination of factor levels. If pair data is used for estimation it usually is valid also here.
#' @param act Activity definition
#' @param data \code{data.frame} Observed transitions
#' @param statespace0 \code{data.frame} Possible states for each combination of by variables before transition
#' @param statespace1 \code{data.frame} Possible states for each combination of by variables after transition
#' @param factors \code{character} Variables used in the iterative estimation
#' @param by \code{character} Variables not used in the iterative estimation
#' @return Activity definition with statespace
#' @export
build_statespace <- function(act, data, statespace0, statespace1, factors, by) {
  act$factors <- factors
  act$by <- by
  act$statespace <- build_statespace_by(data, statespace0, statespace1, act$statespace0, act$statespace1, factors, by)
  act
}
build_statespace_by <- function(data, statespace0, statespace1, state0, state1, factors, by) {
  if(!all(by %in% names(statespace1))) stop("All by variables should be in statespace")
  if(!all(state0 %in% names(data))) stop("All state variables should be in data")

  if(length(by)==0) return(list(build_statespace1(data, statespace0, statespace1, state0, state1, factors)))

  bydata <- split(data, data[by])
  return(lapply(bydata, function(bydata) {
    bystatespace0 <- merge(bydata[1, by, drop=FALSE], statespace0, by=by, all=FALSE)
    bystatespace1 <- merge(bydata[1, by, drop=FALSE], statespace1, by=by, all=FALSE)
    if(nrow(bystatespace0)==0 || nrow(bystatespace1)==0) {
      cat("No rows in statespace for ")
      print(bydata[by][1,])
      stop("Couldn't figure out statespace.")
    }

    res <- build_statespace1(bydata, bystatespace0, bystatespace1,
                             state0, state1, factors)
    res <- cbind(res, bydata[1, by, drop=FALSE], row.names=NULL)
    res
  }))

}

build_statespace1 <- function(data, statespace0, statespace1, state0, state1, factors) {
  # Assuming that state is a cartesian product of individual variables.
  # Convert everything used by table to factors and save the levels because table only return levels as character.
  varlevels <- NULL
  for(name in names(data)) {
    if(name %in% state1) {
      name1 <- gsub("1$", "", name)
      x <- statespace1[[name1]]
    } else if(name %in% state0) {
      name1 <- gsub("0$", "", name)
      x <- statespace0[[name1]]
    } else {
      x <- data[[name]]
    }
    if(is.factor(data[[name]])) {
      if(any(table(data[[name]])==0))
        warning(paste0("Factor ", name,
                       " is a factor, but not all levels have observations. Consider converting ",
                       name, " to character."))
      varlevels[[name]] <- levels(data[[name]])
    } else {
      varlevels[[name]] <- sort(unique(na.omit(x)))
    }
  }
  if(length(factors)==0) {
    out <- data.frame(varlevels=1)
  } else {
    out <- expand.grid(varlevels[factors], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }
  out$varlevels <- list(varlevels)
  out
}
