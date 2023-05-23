#' Estimate Transition Probabilities from Pairdata
#'
#' Estimate Transition Probabilities from Pairdata
#'
#' Transition probabilities 'move' the forest areas allocated in the cells of state matrix
#' from the initial states in the beginning of a EFDM run step to the end position. This end
#' position will be the initial state of the next EFDM step. Length of a step (=time)
#' in the EFDM run is typically determined by the pairdata. It is the time difference of
#' tree observations. Note that the pairdata can be also constructed from single observation,
#' if the other (pair) observation is estimated or modelled.
#'
#' Each activity needs to have a transition probability. If no pairdata is available,
#' transition probability matrices can be based entirely on a prior defined with
#' expert knowledge.
#'
#' The estimation uses an iterative Bayesian algorithm that is explained in
#' \url{https://github.com/ec-jrc/efdm/blob/master/documents/EFDMinstructions/Seija_Mathematics_behind_EFDM.pdf}.
#' The transition probability estimate is the proportion of observed transitions
#' divided by the number of all transitions from the same starting state.
#' \code{prior} gives the number of prior transitions. For each \code{factor}
#' variable the transitions are counted in classes of all factors before the
#' current factor. The "most important" observations (having all classes right)
#' is counted \code{length(factors)} times, the second most important observations
#' are counted \code{length(factors)-1} times and so on.
#'
#' If pairdata is NULL prior is used by itself.
#'
#' Observations should have 'factor' and 'by' variables and statepairs with 0 and 1
#' suffixes to indicate before and after observations.
#'
#' The estimation algorithm uses information across 'factor' variables, but not
#' across 'by' variables.
#'
#' \code{prior} can either \code{character} or \code{function}.
#' \itemize{
#' \item "nochange" implies that there is one observation where state doesn't change
#' \item "uninformative" when no observations are given all states are as likely
#' \item \code{function(A, dynvar1, dynvar0)} where A is an array of zeros with
#'   dimnames(A) <- c(dynvar1, dynvar0). The function should fill A with the number of
#'   prior transitions and return it.
#' }
#'
#' The statespace is used to fill in the transitions where there are no observations.
#' In special cases the statespace may change as a result of the activity.
#' For example changing the tree species might lead to change in the volume classification used.
#'
#' This function assumes that the dynamic variables are coded as integers starting with 1.
#' Other variables are not restricted.
#'
#' @param dynamicvariables The names of the dynamic variables, \code{character} vector
#' @param pairdata \code{data.frame} Observed transitions
#' @param statespace \code{data.frame} or a \code{list} of two \code{data.frame}s if the state space changes as a result of the activity
#' @param factors \code{character} Variables used by the activity
#' @param by \code{character} Variables that split the state space
#' @param prior function or character
#' @return Transition probabilities as a \code{data.frame}
#' @examples
#'
#' # Estimation can use observed transitions with different levels of factors.
#' statespace <- expand.grid(a=1:2, b=1:2, vol=1:5)
#' pairdata <- data.frame(a=c(1,1,2,2), b=c(1,2,1,2), vol0=c(1,1,1,1), vol1=c(2,3,4,5))
#' state0 <- statespace
#' actprob <- statespace
#' actprob$test <- 1
#' state0$area <- 0
#' state0$area[1] <- 1
#'
#' # With by=c("a", "b") there are two observations: one from prior and the other
#' # from the exact combination of class levels.
#' probs <- estimatetransprobs("vol", pairdata, statespace, by=c("a", "b"), prior="nochange")
#' act1 <- define_activity("test", c("vol"), probs)
#' runEFDM(state0, actprob, list(act1), 1)
#'
#' probs <- estimatetransprobs("vol", pairdata, statespace, factors="a", by="b", prior="nochange")
#' act2 <- define_activity("test", c("vol"), probs)
#' runEFDM(state0, actprob, list(act2), 1)
#'
#' probs <- estimatetransprobs("vol", pairdata, statespace, factors="b", by="a", prior="nochange")
#' act3 <- define_activity("test", c("vol"), probs)
#' runEFDM(state0, actprob, list(act3), 1)
#'
#' # The order of variables in factors argument specifies the order of importance.
#' # Observation that differ in the first variable are counted more times.
#' probs <- estimatetransprobs("vol", pairdata, statespace, factors=c("a", "b"), prior="nochange")
#' act4 <- define_activity("test", c("vol"), probs)
#' runEFDM(state0, actprob, list(act4), 1)
#'
#' probs <- estimatetransprobs("vol", pairdata, statespace, factors=c("b", "a"), prior="nochange")
#' act5 <- define_activity("test", c("vol"), probs)
#' runEFDM(state0, actprob, list(act5), 1)
#'
#'
#' @export
estimatetransprobs <- function(dynamicvariables, pairdata,
                               statespace,
                               factors=character(), by=character(),
                               prior="nochange") {
  act <- define_activity("temp", dynamicvariables)
  if (is.data.frame(statespace)) {
    statespace0 <- statespace1 <- statespace
  } else {
    statespace0 <- statespace[[1]]
    statespace1 <- statespace[[2]]
  }
  act <- build_complex_statespace(act, statespace0, statespace1, factors, by)
  dynvar0 <- act$dynamicvariables0
  dynvar1 <- act$dynamicvariables1
  statespace <- act$statespace

  if(!is.null(pairdata)) {
    requirednames <- c(dynvar0, dynvar1, factors, by)
    missingnames <- setdiff(requirednames, names(pairdata))
    if(length(missingnames) > 0) stop(paste0("Variable(s) '", list(missingnames), "' not present in pairdata."))

    # Check pairdata against statespace
    a <- act$statespace0
    a$check_missing_combinations <- 1
    a <- merge(pairdata, a, by.x = c(factors, by, act$dynamicvariables0),
          by.y = c(factors, by, gsub("0$", "", act$dynamicvariables0)), all.x=TRUE)
    if(anyNA(a$check_missing_combinations)) stop("Pairdata has transitions starting from states that are not in statespace.")

    a <- act$statespace1
    a$check_missing_combinations <- 1
    a <- merge(pairdata, a, by.x = c(factors, by, act$dynamicvariables1),
               by.y = c(factors, by, gsub("1$", "", act$dynamicvariables1)), all.x=TRUE)
    if(anyNA(a$check_missing_combinations)) {
      print("Pairdata has transitions resulting in states that are not in statespace.")
      print(utils::head(a[is.na(a$check_missing_combinations),,drop=FALSE]))
      stop("Pairdata has transitions resulting in states that are not in statespace.")
    }
  }

  pairdata <- pairdata[c(dynvar1, dynvar0, factors, by)]
  processed_rows <- 0

  A <- do.call(rbind, lapply(statespace, function(statespacepart) {
    bydata <- NULL
    priorpart <- prior
    if(!is.null(pairdata)) {
      mergeby <- c(factors, by) # Why not just 'by' here?
      if(length(mergeby)) {
        select <- unique(statespacepart$statespace0[mergeby])
        bydata <- as.data.frame(merge(as.data.table(pairdata), as.data.table(select), all=FALSE, by=mergeby))
        if(is.data.frame(prior)) {
          priorpart <- merge(prior, select, all=FALSE, by=mergeby)[c(act$dynamicvariables0, act$dynamicvariables1)]
        }
      } else {
        bydata <- pairdata
      }
      processed_rows <<- processed_rows + nrow(bydata)
    }
    estimatetransprobs1(bydata, statespacepart, priorpart, dynvar0, dynvar1, factors)
  }))
  if(!is.null(pairdata)) {
    if(processed_rows > nrow(pairdata)) stop("Internal error. Processed pairdata multiple times.")
    if(processed_rows < nrow(pairdata)) warning("Not all pairdata was processed. Was this intentional?")
  }
  A$N <- NULL
  A$nobs <- NULL
  A
}

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
    bys0 <- as.data.frame(merge(as.data.table(byrows[i, , drop=FALSE]), as.data.table(statespace0), by=by, all=FALSE))
    bys1 <- as.data.frame(merge(as.data.table(byrows[i, , drop=FALSE]), as.data.table(statespace1), by=by, all=FALSE))
    list(statespace0=bys0, statespace1=bys1)
  }))
}


# Matrix form of prior should only work with numeric states
# => no need to worry about factors
apply_prior <- function(statespace, prior, dynvar0, dynvar1) {
  v0 <- sapply(dynvar0, function(v) {list(sort(unique(statespace$statespace0[[gsub("0$", "", v)]])))})
  v1 <- sapply(dynvar1, function(v) {list(sort(unique(statespace$statespace1[[gsub("1$", "", v)]])))})

  if(is.function(prior)) {
    A <- array(0, dim = sapply(c(v1, v0), length))
    dimnames(A) <- c(v1, v0)
    prior <- prior(A, v1, v0)
    if(!identical(dim(A), dim(prior))) stop("Prior function didn't return array of right size.")
    if(!identical(dimnames(A), dimnames(prior))) stop("Prior function didn't return array with right dimnames.")
  } else if(prior == "nochange") {
    if(!identical(unname(v0), unname(v1))) stop("nochange prior is only available if statespaces are identical")
    prior <- diag(prod(sapply(v0, length)))
    dim(prior) <- sapply(c(v1, v0), length)
    dimnames(prior) <- c(v1, v0)
  } else if(prior == "uninformative") {
    prior <- array(1e-15, dim = sapply(c(v1, v0), length))
    dimnames(prior) <- c(v1, v0)
  } else {
    stop("Only nochange and uninformative prior is implemented")
  }
  check <- prior
  dim(check) <- c(prod(sapply(v1, length)), prod(sapply(v0, length)))
  if(any(colSums(check)==0)) stop("Prior gave zero transition probabilities when starting from some states. Please check your prior.")

  out <- expand.grid(c(v1, v0), stringsAsFactors = FALSE)
  out$nobs <- c(prior)
  out <- out[out$nobs != 0,]
  out <- merge(out, statespace$statespace0, by.x=dynvar0, by.y=gsub("0$", "", dynvar0), sort=FALSE)
  out
}

estimatetransprobs1 <- function(pairdata, statespace, prior, dynvar0, dynvar1, factors) {

  nobs <- if(is.data.frame(prior)) {
    prior
  } else {
    apply_prior(statespace, prior, dynvar0, dynvar1)
  }

  if(!is.null(pairdata) && nrow(pairdata) > 0) {
    # Reorder columns for iterativebayesian estimation
    usednames <- c(dynvar1, dynvar0, factors)

    pairdatafact <- pairdata[usednames]
    # Factors need all levels counted to borrow information across levels
    # dynamic variables only those with observations.
    # factors need to be converted to factor, because different 'by' classes can have different factor levels.
    for(n in factors) {
      pairdatafact[[n]] <- factor(pairdata[[n]], levels=sort(unique(statespace$statespace0[[n]])))
    }
    nobs2 <- iterativebayesestimation(pairdatafact, 0, length(dynvar0), length(dynvar1))
    # do as.data.frame(nobs2), but without converting numberic variables to factors
    vars <- c(sapply(c(dynvar1, dynvar0), function(x) list(sort(unique(pairdata[[x]])))),
      sapply(factors, function(x) list(sort(unique(statespace$statespace0[[x]])))))
    grid <- expand.grid(vars)#, stringsAsFactors = FALSE)
    grid$nobs <- c(nobs2)
    grid <- grid[grid$nobs != 0,]
    # Add by variables
    for(n in setdiff(names(pairdata), usednames)) {
      grid[[n]] <- pairdata[[n]][1]
    }
    if(length(factors)) {
      grid <- merge(grid, unique(statespace$statespace0[factors]), by=factors, all=FALSE)
    }
    nobs <- merge(nobs, grid, by=names(pairdata), suffixes=c("", ".y"), all=TRUE)
    na20 <- function(x) {x[is.na(x)] <- 0; x}
    nobs$nobs <- na20(nobs$nobs) + na20(nobs$nobs.y)
    nobs$nobs.y <- NULL
    if(length(grep("\\.y$", names(nobs)))) stop(paste0("Internal error. Extra variables after transition matrix estimation."))
  }

  N <- aggregate(list(N=nobs$nobs), nobs[c(dynvar0, factors)], sum)
  p <- merge(nobs, N, by=c(dynvar0, factors))
  p$prob <- p$nobs / p$N
  sump <- nrow(statespace$statespace0)
  if(!isTRUE(all.equal(sum(p$prob),sump))) {
    stop(paste0("Sum of estimated transition probabilities is ", sum(p$prob), " but should be ", sump))
  }
  return(p)
}

iterativebayesestimation <- function(pairdata, prior, ls0, ls1) {
  n <- prior
  # Count number of observations with different factors
  # The observations with more specific factor combinations are counted more
  # times than those with less specific combinations.
  for (i in (ls0+ls1):length(pairdata))
  {
    tmp <- table(pairdata[1:i])
    # Since the new dimension is the last dimension it is possible to use
    # recycling to add the previous array to each level of the new dimension
    tmp[] <- c(tmp) + c(n)
    n <- tmp
  }
  n
}
