

do_activity <- function(state, act) {
  requirednames <- gsub("0$", "", act$statespace0)
  missingnames <- setdiff(requirednames, names(state))
  if(length(missingnames) > 0) stop("Variable '", list(missingnames),
                                    "' required by activity '", act$actname, "' not present in state.")


  A <- extract_transitions(act)
  names1 <- setdiff(names(A), c(act$statespace1, "prob"))
  res <- merge(state, A, by.x = gsub("0$", "", names1), by.y=names1, all.x=FALSE, all.y=FALSE)
  res$area <- res$area * res$prob
  res$prob <- NULL
  resbefore <- res
  for(n in act$statespace0)
    res[[gsub("0$", "", n)]] <- NULL
  for(n in act$statespace1) {
    i <- match(n, names(res))
    names(res)[i] <- gsub("1$", "", names(res)[i])
    resbefore[[n]] <- NULL
  }
  resbefore <- aggregate(resbefore["area"], resbefore[names(resbefore) != "area"], sum)
  res <- aggregate(res["area"], res[names(res) != "area"], sum)

  if(!isTRUE(all.equal(sum(resbefore$area), sum(res$area)))) {
    stop(paste0("Some area (", sum(resbefore$area), "!=", sum(res$area), ") was lost in the computation for act ", act$actname))
  }
  list(before=resbefore, after=res)
}

check_activities <- function(acts) {
  actnames <- sapply(acts, function(x) x$actprobname)
  splitacts <- split(acts, actnames)
  if(any(sapply(splitacts, check_activity_split))) stop("Duplicated activies found.")
}

#' Run Markov Forest Dynamics Model
#'
#'
#'
#'
#' @param state0 \code{data.frame} Initial state
#' @param actprob \code{data.frame} Activation probabilities
#' @param activities \code{list} A list of activities
#' @param n \code{integer} Number of time steps required
#' @return \code{data.frame} State of each time step
#' @importFrom stats aggregate na.omit
#' @importFrom utils head
#' @export
runEFDM <- function(state0, actprob, activities, n) {
  totalarea <- sum(state0$area)
  state <- state0
  factornames <- setdiff(names(state), "area")
  actnames <- setdiff(names(actprob), names(state))

  if(anyDuplicated(actprob[factornames]))
    stop("Duplicated rows in actprob.")

  if(any(actprob[actnames] < 0)) stop("Probabilities should be positive.")

  maxdiff <- max(abs(rowSums(actprob[actnames])-1))
  if(maxdiff > 1e-15) stop(paste0("Not all activity probabilities sum to 1. Maximum absolute difference ", maxdiff))

  check_activities(activities)
  beforeactivity <- NULL
  for(i in 0:n) {
    m <- merge(state, actprob, by=factornames, all.x=TRUE)
    if(any(is.na(m[actnames]))) {
      print("No activation probability for:")
      print(head(m[is.na(m[length(m)]),]))
      print("...")
      stop("Stopping")
    }
    res <- list()
    for(act in activities) {
      m2 <- m
      m2$area <- m2$area*m2[[act$actprobname]]
      m2 <- m2[c(factornames, "area")]
      m2 <- m2[m2$area != 0,, drop=FALSE]
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

    allres <- do.call(rbind, res)
    state <- aggregate(allres["area"], allres[factornames], sum)
    newtotalarea <- sum(state$area)
    if(!isTRUE(all.equal(totalarea, newtotalarea)))
      warning(paste("Starting with", totalarea, "area ended up with", newtotalarea, " area."))
  }
  do.call(rbind, beforeactivity)
}
