#' Priors for \code{\link{estimatetransprobs}}
#'
#' Priors for \code{\link{estimatetransprobs}}
#'
#' \code{prior_ff} moves the forest area to the smallest classes of the given
#' dynamic variables of the forest stratum.
#'
#' \code{prior_grow} moves the forest to another class given by increasing
#' \code{variable} by \code{howmuch}.
#'
#' @param variable Name of the variable to grow
#' @param howmuch Amount of growth
#'
#' @return Return value is used by \code{\link{estimatetransprobs}} to provide
#' prior information on the transition probabilities.
#'
#' @examples
#' statespace <- expand.grid(a=1:2, b=1:2, vol=1:15, age=1:35)
#' act <- define_activity("test", c("vol", "age"))
#' act <- build_statespace(act, statespace, by=c("a", "b"))
#' act1 <- estimatetransprobs(act, NULL, prior_ff())
#' act2 <- estimatetransprobs(act, NULL, prior_grow("age"))
#' @export
prior_ff <- function() {
  function(A, state1, state0) {
    s0 <- expand.grid(sapply(state0, function(x) list(seq_along(x))))

    for(n in names(state1)) {
      s0[[n]] <- 1
    }
    A[as.matrix(s0[c(names(state1), names(state0))])] <- 1
    A
  }
}

#' @rdname prior_ff
#' @export
prior_grow <- function(variable, howmuch=1) {
  function(A, state1, state0) {
    stopifnot(identical(unname(state1), unname(state0)))
    vi <- match(paste0(variable, "1"), names(state1))
    if(is.na(vi)) stop(paste0("Variable '", variable, "' not in statespace. Invalid prior specification."))
    v <- state1[[vi]]
    l <- length(v)
    grid0 <- expand.grid(state1)
    state1[[vi]] <- state1[[vi]][pmax(1, pmin((1:l)+howmuch, l))]
    grid1 <- expand.grid(state1)
    A[as.matrix(cbind(grid1, grid0))] <- 1
    A
  }
}
