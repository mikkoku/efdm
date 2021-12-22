#' Define an activity
#'
#' Define an activity
#'
#' The set of activities in EFDM defines all possible alternatives for a forest
#' stratum to develop during a scenario run step. Therefore activities are not only
#' limited to forest treatments and management actions such as thinnings and final fellings
#' but should also include 'no management' i.e. growth, if applicable. An activity may
#' also be something else affecting the development, for example, a forest hazard:
#' snow, wind, drought, pest damage etc.
#'
#' \code{name} is used as the name of the activity in the result of \code{\link{runEFDM}}
#' and it is used to identify the activity probability.
#' If multiple activities should be reported under a common name, \code{actprobname} can
#' be used to specify the activity name in activity probability data.
#'
#' \code{dynamicvariables} are those (state space) variables which change as
#' a result of the activity
#' Typically an activity affects on the age, volume or stem count of the
#' forest, but an activity may also, for example, change land-use and then
#' a variable related to land-use categories is essential.
#'
#' The \code{transprobs} should be a \code{data.frame} with columns:
#' \itemize{
#' \item pairs of dynamic variables. For example, if the dynamic variables are
#'   \code{vol} and \code{age}, then transprobs should have columns
#'   \code{age0}, \code{age1}, \code{vol0}, \code{vol1}.
#' \item the probability of transition \code{prob}
#' \item and possibly other variables affecting the transition probabilities.
#'   For example \code{soiltype}, \code{treespecies}, etc.
#' }
#' \code{\link{estimatetransprobs}} can be used to estimate transition probabilities
#' from a similar set of pair data of observations.
#'
#' @return An activity definition
#'
#' @param name Name of activity used in reporting
#' @param dynamicvariables Names of variables where changes happen
#' @param transprobs Transition probabilities (\code{data.frame})
#' @param probname (optional) Name of activity in activity probabilities data
#' @examples
#' define_activity("nomanagement", c("vol", "age"))
#' @export
define_activity <- function(name, dynamicvariables, transprobs=NULL, probname=name) {
  act <- list(actname=name,
       actprobname=probname,
       dynamicvariables0=paste0(dynamicvariables, "0", recycle0=TRUE),
       dynamicvariables1=paste0(dynamicvariables, "1", recycle0=TRUE))
  if (!is.null(transprobs))
    efdm::transprobs(act) <- transprobs
  act
}


#' @rdname transprobs
#' @export
`transprobs<-` <- function(act, value) {
  probs <- value
  if(is.null(act$dynamicvariables0) || is.null(act$dynamicvariables1) || is.null(act$actname)) stop("Not a valid activity.")
  if(length(act$dynamicvariables0) == 0)
    stop("Activity without dynamic variables may not have transition probabilities.")
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
#' The \code{value} should be a \code{data.frame} defining the transition probabilities
#' between dynamic variables of the activity. See \code{\link{define_activity}}
#' for details
#' @param act Activity definition
#' @param value \code{data.frame} of transition probabilities
#' @return \code{data.frame} where prob is the transition probability from current
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
