extract_A <- function(act) {
  A <- act$A
  A$nobs <- NULL
  A$N <- NULL
  A
}

aggregate_area <- function(res, by=names(res)) {
  ns <- setdiff(by, "area")
  r <- as.data.table(res)
  # Using area instead of .SD$area appears to be a lot faster
  area <- stop # This is to trick R Check into not warning about the variable "area"
  as.data.frame(r[,list(area=sum(area)), by=ns])
}

#' @importFrom data.table as.data.table .SD
do_activity <- function(state, act) {
  if(length(act$dynamicvariables0) == 0)
    return(list(before=state, after=state))

  A <- extract_A(act)
  requirednames <- c(setdiff(names(A), c("prob", act$dynamicvariables0, act$dynamicvariables1)),
    gsub("0$", "", act$dynamicvariables0))
  missingnames <- setdiff(requirednames, names(state))
  if(length(missingnames) > 0) stop("Variable '", list(missingnames),
                                    "' required by activity '", act$actname, "' not present in state.")


  names1 <- intersect(names(A), names(state))
                        #3c(act$dynamicvariables1, "prob"))
  by.x = c(names1, gsub("0$", "", act$dynamicvariables0))
  by.y = c(names1, act$dynamicvariables0)
  res <- merge(as.data.table(state), as.data.table(A), by.x=by.x, by.y=by.y,
               all.x=FALSE, all.y=FALSE, allow.cartesian=TRUE)
  res$area <- res$area * res$prob
  res$prob <- NULL
  resbefore <- res
  for(n in act$dynamicvariables0)
    res[[gsub("0$", "", n)]] <- NULL
  for(n in act$dynamicvariables1) {
    i <- match(n, names(res))
    names(res)[i] <- gsub("1$", "", names(res)[i])
    resbefore[[n]] <- NULL
  }
  resbefore <- aggregate_area(resbefore)
  res <- aggregate_area(res)

  if(!isTRUE(all.equal(sum(resbefore$area), sum(res$area)))) {
    stop(paste0("Some area (", sum(resbefore$area), "!=", sum(res$area), ") was lost in the computation for act ", act$actname))
  }
  list(before=resbefore, after=res)
}

check_activities <- function(acts, actprob, actnames, actfactors) {
  actnames1 <- sapply(acts, function(x) x$actprobname)
  if(!all(actnames1 %in% actnames)) stop(paste0("Activity probabilities for ", setdiff(actnames1, actnames), " not given."))
  if(!all(actnames %in% actnames1)) stop(paste0("Activity probabilities for ", setdiff(actnames, actnames1), " given, but activities not."))

  P <- actprob
  P$check___ <- 0
  if(anyNA(P)) stop("Actprob should not have NAs.")

  if(any(P[actnames] < 0)) stop("Activity probabilities should be positive.")

  maxdiff <- max(abs(rowSums(P[actnames])-1))
  if(maxdiff > 1e-15) stop(paste0("Not all activity probabilities sum to 1. Maximum absolute difference ", maxdiff))


  # Check states before transition
  for(act in acts) {
    A <- extract_A(act)
    for(v in act$dynamicvariables1)
      A[[v]] <- NULL
    for(v0 in act$dynamicvariables0) {
      v <- gsub("0$", "", v0)
      A[[v]] <- A[[v0]]
      A[[v0]] <- NULL
    }

    A$prob <- NULL
    A <- A[!duplicated(A),,drop=FALSE]
    A$check_new__ <- 1
    P <- merge(P, A, all.x=TRUE, all.y=TRUE)
    P$check___ <- P$check___ + ifelse(is.na(P$check_new__), 0, P[[act$actname]])
    P$check_new__ <- NULL
  }
  check <- is.na(P$check___)
  if(any(check)) {
    print("Activity transitions from state not in actprob:")
    print(utils::head(P[check,,drop=FALSE]))
    stop("Activity transitions from state not in actprob")
  }
  check <- P$check___ == 0
  if(any(check)) {
    print("States with no activity:")
    print(utils::head(P[check,,drop=FALSE]))
    stop("States with no activity")
  }
  P$check_1__ <- P$check___ - 1
  check <- abs(P$check_1__) > 1e-14
  if(any(check)) {
    print("States with multiple activities:")
    print(utils::head(P[check,,drop=FALSE]))
    stop("States with multiple activities")
  }

  # Check states after transition
  for(act in acts) {
    if(length(act$dynamicvariables0) > 0) {
      A <- extract_A(act)
      for(v in act$dynamicvariables0)
        A[[v]] <- NULL
      for(v1 in act$dynamicvariables1) {
        v <- gsub("1$", "", v1)
        A[[v]] <- A[[v1]]
        A[[v1]] <- NULL
      }
      A$prob <- NULL
      A <- A[!duplicated(A),,drop=FALSE]
      P <- merge(P, A, all.x=TRUE, all.y=TRUE)
    }
  }
  check <- is.na(P$check___)
  if(any(check)) {
    print("Activity transitions to state not in actprob:")
    print(utils::head(P[check,,drop=FALSE]))
    stop("Activity transitions to state not in actprob")
  }
}



#' Run European Forestry Dynamics Model
#'
#' Run European Forestry Dynamics Model
#'
#' This is the actual scenario running function, which projects the initial
#' forest state n time steps to the future.
#'
#' An activity is defined by activity name, names of dynamic variables and transition probabilities.
#'
#'
#' @param state0 \code{data.frame} Initial state
#' @param actprob \code{data.frame} Activity probabilities
#' @param activities \code{list} A list of activities
#' @param n \code{integer} Number of time steps required
#' @param check Check input arguments for consistency.
#' @return \code{data.frame} State of each time step divided by activities
#' @importFrom stats aggregate
#' @importFrom utils head
#' @export
runEFDM <- function(state0, actprob, activities, n, check=TRUE) {
  totalarea <- sum(state0$area)
  state <- state0
  # state may include for example plot id which should not be anywhere else
  actfactors <- intersect(names(state), names(actprob))
  actnames <- setdiff(names(actprob), names(state))

  if(length(actfactors)==0) stop("No state variables in actprob.")

  if(check) {
    if(anyDuplicated(actprob[setdiff(names(actprob), actnames)]))
      stop("Duplicated rows in actprob.")

    check_activities(activities, actprob, actnames, actfactors)
  }

  beforeactivity <- NULL
  for(i in 0:n) {
    m <- as.data.frame(merge(as.data.table(state), as.data.table(actprob), by=actfactors, all.x=TRUE, all.y=FALSE, allow.cartesian=TRUE))

    if(any(is.na(m[actnames]))) {
      # This shouldn't happen unless check=FALSE
      print("No activation probability for:")
      print(head(m[is.na(m[length(m)]),]))
      print("...")
      stop("Stopping")
    }
    checkareain <- list()
    checkareaout <- list()
    for(actname in actnames) {
      checkareaout[[actname]] <- 0
    }

    res <- list()
    for(act in activities) {
      m2 <- m
      m2$area <- m2$area*m2[[act$actprobname]]
      m2 <- m2[names(state)]
      m2 <- m2[m2$area != 0,, drop=FALSE]
      checkareain[[act$actprobname]] <- sum(m2$area)
      # m2 can include some area that is not meant for this activity due to
      # activities having different statespace but shared activation probability
      if(nrow(m2) != 0) {
        res2 <- do_activity(m2, act)
        res1 <- res2$after
        if(nrow(res1) != 0)  {
          extranames <- setdiff(names(state), c(names(res1), "area"))
          if(length(extranames) != 0) {
            stop(paste0("Variable '", list(extranames), "' in state before activity '", act$actname, "' but not after."))
          }
          res1$activity <- act$actname
          res1$time <- i
          res <- c(res, list(res1))
          checkareaout[[act$actprobname]] <- checkareaout[[act$actprobname]] + sum(res1$area)
        }
        res1 <- res2$before
        res1 <- res1[!is.na(res1$area),,drop=FALSE]
        if(nrow(res1) != 0)  {
          res1$activity <- act$actname
          res1$time <- i
          beforeactivity <- c(beforeactivity, list(res1))
        }
      }
    }
    # These checks should not be needed because they are made in check_activities
    for(actname in actnames) {
      if(!isTRUE(all.equal(checkareain[[actname]], checkareaout[[actname]])))
        stop(paste0("Activity '", actname, "' was allocated ", checkareain[[actname]], " area, but it processed ", checkareaout[[actname]], "."))
    }

    if(length(res) == 0) {
      stop(paste0("All area was lost after time ", i, "."))
    }
    allres <- do.call(rbind, res)
    state <- aggregate_area(allres, names(state))
    newtotalarea <- sum(state$area)
    if(!isTRUE(all.equal(totalarea, newtotalarea)))
      warning(paste("Starting with", totalarea, "area ended up with", newtotalarea, " area."))
  }
  do.call(rbind, beforeactivity)
}
