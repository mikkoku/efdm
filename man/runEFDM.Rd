% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.r
\name{runEFDM}
\alias{runEFDM}
\title{Run European Forestry Dynamics Model}
\usage{
runEFDM(state0, actprob, activities, n, check = TRUE)
}
\arguments{
\item{state0}{\code{data.frame} Initial state}

\item{actprob}{\code{data.frame} Activity probabilities}

\item{activities}{\code{list} A list of activities}

\item{n}{\code{integer} Number of time steps required}

\item{check}{Check input arguments for consistency.}
}
\value{
\code{data.frame} State of each time step divided by activities
}
\description{
Run European Forestry Dynamics Model
}
\details{
This is the actual scenario running function, which projects the initial
forest state n time steps to the future.

Activities are defined using \code{\link{define_activity}}.
}
