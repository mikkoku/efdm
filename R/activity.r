#' Define activity
#'
#' @param name Name of activity used in reporting
#' @param dynamicvariables Names of variables where changes happen
#' @param probname (optional) Name of activity in activity probabilities data
#' @examples
#' define_activity("nomanagement", c("vol", "age"))
#' @export
define_activity <- function(name, dynamicvariables, probname=name) {
  list(actname=name,
       actprobname=probname,
       dynamicvariables0=paste0(dynamicvariables, "0"),
       dynamicvariables1=paste0(dynamicvariables, "1"))
}

#' Add a statespace to an activity
#'
#'
#' Statespace is the collection of strata.
#'
#' \code{factors} and \code{by} variables are used by \code{\link{estimatetransprobs}}
#' to estimate transition probabilities for missing strata. Observations with
#' different levels of \code{factors} variables are used where as observations
#' with different levels of \code{by} variables are never used.
#'
#' If statespace is changing as a result of the activity see \code{\link{build_complex_statespace}}.
#' @param act Activity definition
#' @param statespace \code{data.frame}
#' @param factors \code{character} Variables used by the activity
#' @param by \code{character} Variables that split the statespace
#' @return Activity definition with statespace
#' @examples
#'
#' statespacepine <- expand.grid(species="pine", vol=1:5, age=1:3, stringsAsFactors=FALSE)
#' statespacespruce <- expand.grid(species="spruce", vol=1:4, age=1:3, stringsAsFactors=FALSE)
#' statespace <- rbind(statespacepine, statespacespruce)
#' act <- define_activity("nomanagement", c("vol", "age"))
#' act <- build_statespace(act, statespace, by="species")
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
#' if(sum(sapply(act$statespace, function(x) nrow(x$statespace0)))!=length(ii)) {
#'   stop(paste("build_statespace failed", list(ii)))
#' }
#' }
#' @export
build_statespace <- function(act, statespace, factors=character(), by=character()) {
  # rename factors => class
  build_complex_statespace(act, statespace, statespace, factors, by)
}
#' Add a statespace to an activity
#'
#'
#' In a complex statespace it is possible to change statespaces with an activity.
#' @inheritParams build_statespace
#' @param act Activity definition
# @param data \code{data.frame} Observed transitions
#' @param statespace0 \code{data.frame} Statespace before transition
#' @param statespace1 \code{data.frame} Statespace after transition
#' @seealso build_state_space
#' @return Activity definition with statespace
#' @export
build_complex_statespace <- function(act, statespace0, statespace1, factors=character(), by=character()) {
  extranames <- setdiff(names(statespace0), c(factors, by, gsub("0$", "", act$dynamicvariables0)))
  if(length(extranames)) stop(paste0("Variables '", list(extranames), "' in statespace not used by activity."))
  act$factors <- factors
  act$by <- by
  act$statespace0 <- statespace0
  act$statespace1 <- statespace1
  act$statespace <- build_statespace_by(statespace0, statespace1, act$dynamicvariables0, act$dynamicvariables1, factors, by)
  act
}
build_statespace_by <- function(statespace0, statespace1, dynamicvariables0, dynamicvariables1, factors, by) {
  if(!all(by %in% names(statespace1))) stop("All by variables should be in statespace1")
  if(!all(gsub("0$", "", dynamicvariables0) %in% names(statespace0))) stop("All dynamic variables should be in statespace0")
  if(!all(gsub("1$", "", dynamicvariables1) %in% names(statespace1))) stop("All dynamic variables should be in statespace1")
  if(!all(by %in% names(statespace0))) stop("All by variables should be in statespace0")

  if(length(by)==0) return(list(list(statespace0=statespace0, statespace1=statespace1)))
  byrows <- unique(rbind(statespace0[by], statespace1[by]))
  return(lapply(1:nrow(byrows), function(i) {
    #TODO: Change merge to a simple subset
    # merge has an interesting view on "sort"ing
    bys0 <- merge(byrows[i, , drop=FALSE], statespace0, by=by, all=FALSE)
    bys1 <- merge(byrows[i, , drop=FALSE], statespace1, by=by, all=FALSE)
    list(statespace0=bys0, statespace1=bys1)
  }))
}

#' Extract transition probabilities from an activity
#'
#' @param act Activity
#' @return \code{data.frame} where prob=nobs/N is the transition probability from current
#' state (with suffix 0) to next state (with suffix 1).
#' @export
extract_transitions <- function(act) {
  act$A
}
