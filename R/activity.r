#' Define an activity
#'
#' Define an activity
#'
#' The set of activies in EFDM defines all possible alternatives for a forest
#' stratum to develop during a scenario run step. Therefore activities are not only
#' limited to forest treatments and management actions such as thinnings and final fellings
#' but should also include 'no management' i.e. growth, if applicable. An activity may
#' also be something else affecting the development, for example, a forest hazard:
#' snow, wind, drought, pest damage etc.
#'
#' An activity is defined with this function. A name, which is henceforth
#' used in the EFDM R project when refering to the activity, is given. In addition,
#' the (stetaspace) variables which are affected by the activity are named.
#' Typically an activity affects on the age, volume or stem count of the
#' forest, but an activity may also, for example, change land-use and then
#' a variable related to land-use categories is essential. If the activity name does
#' not match the name in the activaty probability data set, those can be linked here.
#'
#' Defining an activity is the first step, which will be followed by
#' \itemize{
#' \item \code{\link{build_statespace}} or \code{\link{build_complex_statespace}} and
#'   \code{\link{estimatetransprobs}} if using pairdata
#' \item \code{\link{transprobs<-}} if not using pairdata
#' }
#' until an activity is fully applicable in runEFDM.
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
#' Add a statespace to an activity
#'
#' Statespace is the collection of strata. When it is added to an activity defined by define_activity
#' the following \code{\link{estimatetransprobs}} function has sufficient information to estimate transition
#' probabilties.
#'
#' \code{factors} and \code{by} variables are used by \code{\link{estimatetransprobs}}
#' to estimate transition probabilities for strata. Observations with
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
#' @export
build_statespace <- function(act, statespace, factors=character(), by=character()) {
  # rename factors => class
  build_complex_statespace(act, statespace, statespace, factors, by)
}
#' Add a statespace to an activity
#'
#' Add a statespace to an activity
#'
#' In a complex statespace it is possible to change statespaces with an activity.
#' Since statespace is the collection of classes of variables this means that the
#' classification changes.
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

#' @rdname transprobs
#' @export
`transprobs<-` <- function(act, value) {
  probs <- value
  if(is.null(act$dynamicvariables0) || is.null(act$dynamicvariables1) || is.null(act$actname)) stop("Not a valid activity.")
  if(!("prob" %in% names(probs))) stop("probs should have a column named 'prob'")
  vars <- c(act$dynamicvariables1, act$dynamicvariables0)
  if(!all(vars %in% names(probs)))
    stop(paste0("dynamicvariables ", list(setdiff(vars, names(probs)))), " not found in 'probs'")

  aggregatenames <- setdiff(names(probs), c(act$dynamicvariables1, c("prob", "N", "nobs")))
  totalprobs <- aggregate(probs['prob'], probs[aggregatenames], FUN=sum)

  maxdiff <- max(abs(totalprobs$prob-1))
  if(maxdiff > 1e-15) stop(paste0("Not all probabilities sum to 1. Maximum absolute difference ", maxdiff))

  act$A <- probs
  act
}

#' Transition probabilities of an activity
#'
#' Functions to get or set the transition probabilities of an activity
#'
#' The \code{probs} should contain
#' \itemize{
#' \item dynamic variables in the activity
#' \item the probability of transition \code{prob}
#' \item other variables affecting the transition probabilities.
#' }
#' @param act Activity definition
#' @param value \code{data.frame} of transition probabilities
#' @return \code{data.frame} where prob=nobs/N is the transition probability from current
#' state (with suffix 0) to next state (with suffix 1).
#' @examples
#' act1 <- define_activity("test", c("vol"))
#' transprobs(act1) <- data.frame(vol0 = 1:5, vol1=c(2:5, 5), prob=1)
#' transprobs(act1)
#'
#' if(require("ggplot2")) {
#'   ggplot(transprobs(act1)) + geom_raster(aes(x=vol0, y=vol1, fill=prob))
#' }
#' @export
transprobs <- function(act) {
  act$A
}
#' @rdname transprobs
#' @export
extract_transitions <- function(act) {
  warning("extract_transitions replaced by transprobs.")
  transprobs(act)
}
