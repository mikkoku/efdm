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
#' Statespace is the collection of levels for each state variable.
#'
#'
#' If statespace is changing as a result of the activity see \code{\link{build_complex_statespace}}.
#' @param act Activity definition
#' @param statespace \code{data.frame}
#' @param factors \code{character} Variables used in the iterative estimation
#' @param by \code{character} Variables not used in the iterative estimation
#' @return Activity definition with statespace
#' @examples
#'
#' statespacepine <- expand.grid(species="pine", vol=1:5, age=1:3, stringsAsFactors=FALSE)
#' statespacespruce <- expand.grid(species="spruce", vol=1:4, age=1:3, stringsAsFactors=FALSE)
#' statespace <- rbind(statespacepine, statespacespruce)
#' act <- define_activity("nomanagement", c("vol", "age"))
#' act <- build_statespace(act, statespace)
#'
#' \dontshow{
#' statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3, stringsAsFactors=FALSE)
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"))
#' stopifnot(nrow(act$statespace[[1]])==3)
#'
#' statespace <- expand.grid(a=1:2, b=c("f", "g"), c=c(4,9), vol=1:3, stringsAsFactors=FALSE)
#' ii <- sample(nrow(statespace), sample(nrow(statespace), 1))
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, statespace[ii,], factors=c("a", "b"), by="c")
#' if(sum(sapply(act$statespace, nrow))!=nrow(unique(statespace[ii,c("a", "b", "c")]))) {
#'   stop(paste("build_statespace failed", list(ii)))
#' }
#' }
#' @export
build_statespace <- function(act, statespace, factors=character(), by=character()) {
  statespace0 <- statespace1 <- data <- statespace
  for(name in c(act$statespace0, act$statespace1)) {
    data[[name]] <- data[[substr(name, 1, nchar(name)-1)]]
  }
  build_complex_statespace(act, data, statespace0, statespace1, factors, by)
}
#' Add a statespace to an activity
#'
#'
#' In a complex statespace it is possible to change statespaces with an activity.
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
build_complex_statespace <- function(act, data, statespace0, statespace1, factors=character(), by=character()) {
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
    out <- unique(data[factors])
  }
  out$varlevels <- list(varlevels)
  out
}
