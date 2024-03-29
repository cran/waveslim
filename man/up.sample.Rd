% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/up.sample.R
\name{up.sample}
\alias{up.sample}
\title{Upsampling of a vector}
\usage{
up.sample(x, f, y = NA)
}
\arguments{
\item{x}{vector of observations}

\item{f}{frequency of upsampling; e.g, 2, 4, etc.}

\item{y}{value to upsample with; e.g., NA, 0, etc.}
}
\value{
A vector twice its length.
}
\description{
Upsamples a given vector.
}
\references{
Any basic signal processing text.
}
\author{
B. Whitcher
}
\keyword{ts}
