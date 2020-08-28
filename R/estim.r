#' Estimate Transition Probabilities from Pairdata
#'
#'
#' The estimation uses an iterative Bayesian algorithm that is explained in
#' \url{https://github.com/ec-jrc/efdm/blob/master/documents/EFDMinstructions/Seija_Mathematics_behind_EFDM.pdf}.
#'
#' If pairdata is NULL prior is used by itself.
#'
#' @param act Activity definition with statespace
#' @param pairdata \code{data.frame} Observed transitions
#' @param prior function or character
#' @return Activity definition with transition probabilities
#' @export
estimatetransprobs <- function(act, pairdata, prior) {
  stopifnot(!is.null(act$statespace))
  state0 <- act$statespace0
  state1 <- act$statespace1
  factors <- act$factors
  by <- act$by
  statespace <- act$statespace

  act$transmat <- do.call(rbind, lapply(statespace, function(statespacepart) {
    bydata <- NULL
    if(!is.null(pairdata))
      bydata <- merge(pairdata, statespacepart[by], all=FALSE, by=by)
    estimatetransprobs1(bydata, statespacepart, prior, state0, state1, factors)
  }))

  act
}

apply_prior <- function(statespace, prior, state0, state1) {
  varlevels <- statespace[["varlevels"]][[1]]
  v0 <- varlevels[state0]
  v1 <- varlevels[state1]

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
  prior
}

#' @importFrom stats na.omit
estimatetransprobs1 <- function(pairdata, statespace, prior, state0, state1, factors) {
  usednames <- c(state1, state0, factors)
  # Check that other variables are constant
  for(name in setdiff(names(pairdata), usednames)) {
    if(!all(pairdata[[name]] == pairdata[[name]][1]))
      stop(paste("Don't know what to do with '", name, "' variable."))
  }

  varlevels <- statespace[["varlevels"]][[1]]

  s <- apply_prior(statespace, prior, state0, state1)

  if(!is.null(pairdata)) {
    pairdata <- pairdata[,usednames, drop=FALSE]

    for(name in names(pairdata)) {
      if(!is.factor(pairdata[[name]])) {
        pairdata[[name]] <- factor(pairdata[[name]], levels=varlevels[[name]])
      }
    }


    s <- iterativebayesestimation(pairdata, s, length(state0), length(state1))
  }

  # Convert s to array and add names to dim
  class(s) <- "array"
  dims <- dimnames(s)
  names(dim(s)) <- names(dimnames(s))
  dimnames(s) <- dims
  transmat <- s

  out <- statespace

  if(length(factors) == 0) {
    out$transmat[[1]] <- transmat
  } else {
    factordims <- dimnames(transmat)[factors]
    grid <- expand.grid(factordims, KEEP.OUT.ATTRS = FALSE)
    out$transmat <- list(NULL)
    for(j in 1:nrow(out)) {
      d <- c(dims[state1], dims[state0], grid[j,])[seq_along(dim(s))]
      # Drop singleton dimensions only for "grid" variables
      A <- do.call("[", c(list(transmat), d, list(drop=FALSE)))
      dn <- dimnames(A)
      dim(A) <- dim(A)[1:length(c(state0, state1))]
      dimnames(A) <- dn[1:length(dim(A))]
      out$transmat[[j]] <- A
    }
  }
  out
}

iterativebayesestimation <- function(pairdata, prior, ls0, ls1) {
  s <- prior
  for (i in (ls0+ls1):length(pairdata))
  {
    tmp <- table(pairdata[1:i])
    # Since the new dimension is the last dimension it is possible to use
    # recycling to add the previous array to each level of the new dimension
    tmp[] <- c(tmp) + c(s)
    s <- tmp
  }

  # s contains all the numerator values
  # to obtain the denominator values the dimensions corresponding to state1
  # are summed over.
  N <- apply(s, (ls1+1):length(dim(s)), sum)
  sweep(s, (ls1+1):length(dim(s)), N, '/')
}
