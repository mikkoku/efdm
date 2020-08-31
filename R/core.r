

do_activity <- function(state, act) {
  res <- NULL
  resbefore <- NULL

  grid <- act$transmat
  for(i in 1:nrow(grid)) {
    currentfactors <- grid[i,,drop=FALSE]
    row.names(currentfactors) <- NULL
    transmat <- currentfactors$transmat[[1]]
    varlevels <- currentfactors$varlevels[[1]]
    currentfactors$transmat <- NULL
    currentfactors$varlevels <- NULL

    statespace0 <- expand.grid(varlevels[act$statespace0], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    statespace1 <- expand.grid(varlevels[act$statespace1], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

    # Select State corresponding to current fixed factors
    m2 <- merge(state, currentfactors, all=FALSE)

    statespace0$seq_used_for_sorting <- 1:nrow(statespace0)

    # Add rows for all combinations of statespace0
    # Remove rows not in statespace0
    m2 <- merge(m2, statespace0, all.x=FALSE, all.y=TRUE, by.x=gsub("0$", "", act$statespace0), by.y=act$statespace0)

    # m2 contains the area that this activity is applied to.
    resbefore <- rbind(resbefore, m2)#[(!is.na(m2$area)) && (m2$area != 0),, drop=FALSE])

    # Restore the ordering in statespace0
    m2 <- m2[order(m2$seq_used_for_sorting),, drop=FALSE]
    m2$seq_used_for_sorting <- NULL

    x <- m2$area
    x[is.na(x)] <- 0
    A <- transmat
    dim(A) <- c(nrow(statespace1), nrow(statespace0))
    y <- A %*% x

    names(statespace1) <- gsub("1$", "", names(statespace1))
    m3 <- cbind(currentfactors, statespace1)
    m3$area <- c(y)
    m3 <- m3[m3$area != 0,, drop=FALSE]

    if(!isTRUE(all.equal(sum(y), sum(x)))) {
      stop(paste0("Some area (", sum(y) - sum(x), ") was lost in the computation for act ", act$actname))
    }
    res <- rbind(res, m3)
  }
  resbefore$seq_used_for_sorting <- NULL
  list(before=resbefore, after=res)
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
  #states <- list() #list(list(sum=state0))
  beforeactivity <- NULL
  afteractivity <- NULL
  for(i in 1:n) {
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
        #beforeactivity <- c(beforeactivity, list(m2))
        res2 <- do_activity(m2, act)
        res1 <- res2$after
        if(nrow(res1) != 0)  {
          if(length(setdiff(names(state), c(names(res1), "area"))) != 0) {
            stop("Missing factor")
          }
          res1$activity <- act$actname
          res1$time <- i
          res <- c(res, list(res1))
          afteractivity <- c(afteractivity, list(res1))
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
    #states[[i]] <- c(list(sum=state), res)
    newtotalarea <- sum(state$area)
    if(!isTRUE(all.equal(totalarea, newtotalarea)))
      warning(paste("Starting with", totalarea, "area ended up with", newtotalarea, " area."))
  }
  #states
  list(before=do.call(rbind, beforeactivity), after=do.call(rbind, afteractivity))
}
