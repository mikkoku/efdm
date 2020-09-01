#' Estimate Transition Probabilities from Pairdata
#'
#'
#' The estimation uses an iterative Bayesian algorithm that is explained in
#' \url{https://github.com/ec-jrc/efdm/blob/master/documents/EFDMinstructions/Seija_Mathematics_behind_EFDM.pdf}.
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
#' \item \code{function(A, state1, state0)} where A is an array of zeros with
#'   dimnames(A) <- c(state1, state0). The function should fill A with observed
#'   transitions and return it.
#' }
#'
#' @param act Activity definition with statespace
#' @param pairdata \code{data.frame} Observed transitions
#' @param prior function or character
#' @return Activity definition with transition probabilities
#' @examples
#' statespace <- expand.grid(ds=c("sp", "pi"), vol=1:3)
#' pairdata <- data.frame(ds=c("sp", "pi"), vol0=c(1,1), vol1=c(2,3))
#' actprob <- state0 <- statespace
#' actprob$test <- 1
#' state0$area <- c(1,1,0,0,0,0)
#'
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, statespace, factors=c("ds"))
#' act1 <- estimatetransprobs(act, pairdata, "nochange")
#' runEFDM(state0, actprob, list(act1), 1)
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, statespace, by=c("ds"))
#' act2 <- estimatetransprobs(act, pairdata, "nochange")
#' runEFDM(state0, actprob, list(act2), 1)
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, statespace, by=c("ds"))
#' growth_model <- function(A, state1, state0) {
#'   A[1,1] <- 0.2; A[2,1] <- 0.7; A[3,1] <- 0.1
#'   A[2,2] <- 0.1; A[3,2] <- 0.9
#'   A[3,3] <- 1
#'   A
#' }
#' act3 <- estimatetransprobs(act, NULL, growth_model)
#' runEFDM(state0, actprob, list(act3), 2)
#'
#' statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3, stringsAsFactors=FALSE)
#' pairdata <- data.frame(ds=c("sp", "pi"), vol0=c(1,1), vol1=c(2,3), region=c("n", "s"))
#' actprob <- state0 <- statespace
#' actprob$test <- 1
#' state0$area <- c(1,1,1,1, 0,0,0,0, 0,0,0,0)
#'
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, subset(statespace, ds=="sp"&region=="n"), factors=c("ds", "region"))
#' act1 <- estimatetransprobs(act, NULL, "nochange")
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"))
#' act2 <- estimatetransprobs(act, NULL, "nochange")
#' runEFDM(state0, actprob, list(act1, act2), 1)
#'
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, subset(statespace, ds=="sp"&region=="n"), factors=c("ds", "region"))
#' act1 <- estimatetransprobs(act, pairdata, "nochange")
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, subset(statespace, !(ds=="sp"&region=="n")), factors=c("ds", "region"))
#' act2 <- estimatetransprobs(act, pairdata, "nochange")
#' runEFDM(state0, actprob, list(act1, act2), 10)
#'
#' statespace <- expand.grid(ds=c("sp", "pi"), region=c("n", "s"), vol=1:3, stringsAsFactors=FALSE)
#' statespace$split <- statespace$ds == "sp" & statespace$region == "n"
#' pairdata <- data.frame(ds=c("sp", "pi"), vol0=c(1,1), vol1=c(2,3), region=c("n", "s"))
#' pairdata$split <- pairdata$ds == "sp" & pairdata$region == "n"
#' actprob <- state0 <- statespace
#' actprob$test <- 1
#' state0$area <- c(1,1,1,1, 0,0,0,0, 0,0,0,0)
#'
#' act <- define_activity("test", c("vol"))
#' act <- build_statespace(act, statespace, factors=c("region"), by=c("split", "ds"))
#' act1 <- estimatetransprobs(act, pairdata, "nochange")
#' runEFDM(state0, actprob, list(act1), 10)
#' @export
estimatetransprobs <- function(act, pairdata, prior) {
  stopifnot(!is.null(act$statespace))
  state0 <- act$statespace0
  state1 <- act$statespace1
  factors <- act$factors
  by <- act$by
  statespace <- act$statespace

  processed_rows <- 0

  act$transmat <- do.call(rbind, lapply(statespace, function(statespacepart) {
    bydata <- NULL
    if(!is.null(pairdata)) {
      bydata <- merge(pairdata, statespacepart[c(factors,by)], all=FALSE, by=c(factors,by))
      processed_rows <<- processed_rows + nrow(bydata)
    }
    estimatetransprobs1(bydata, statespacepart, prior, state0, state1, factors)
  }))
  if(!is.null(pairdata)) {
    if(processed_rows > nrow(pairdata)) stop("Internal error. Processed pairdata multiple times.")
    if(processed_rows < nrow(pairdata)) warning("Not all pairdata was processed. Was this intentional?")
  }
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
  check_unused_variables(pairdata, usednames)

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
