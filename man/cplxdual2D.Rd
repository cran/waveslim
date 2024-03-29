% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cplxdual2D.R
\name{cplxdual2D}
\alias{cplxdual2D}
\alias{icplxdual2D}
\title{Dual-tree Complex 2D Discrete Wavelet Transform}
\usage{
cplxdual2D(x, J, Faf, af)

icplxdual2D(w, J, Fsf, sf)
}
\arguments{
\item{x}{2D array.}

\item{J}{number of stages.}

\item{Faf}{first stage analysis filters for tree i.}

\item{af}{analysis filters for the remaining stages on tree i.}

\item{w}{wavelet coefficients.}

\item{Fsf}{last stage synthesis filters for tree i.}

\item{sf}{synthesis filters for the preceeding stages.}
}
\value{
For the analysis of \code{x}, the output is \item{w}{wavelet
coefficients indexed by \code{[[j]][[i]][[d1]][[d2]]}, where
\eqn{j=1,\ldots,J} (scale), \eqn{i=1} (real part) or \eqn{i=2} (imag part),
\eqn{d1=1,2} and \eqn{d2=1,2,3} (orientations).} For the synthesis of
\code{w}, the output is \item{y}{output signal.}
}
\description{
Dual-tree complex 2D discrete wavelet transform (DWT).
}
\examples{

\dontrun{
## EXAMPLE: cplxdual2D
x = matrix(rnorm(32*32), 32, 32)
J = 5
Faf = FSfarras()$af
Fsf = FSfarras()$sf
af = dualfilt1()$af
sf = dualfilt1()$sf
w = cplxdual2D(x, J, Faf, af)
y = icplxdual2D(w, J, Fsf, sf)
err = x - y
max(abs(err))
}

}
\seealso{
\code{\link{FSfarras}}, \code{\link{farras}}, \code{\link{afb2D}},
\code{\link{sfb2D}}.
}
\author{
Matlab: S. Cai, K. Li and I. Selesnick; R port: B. Whitcher
}
\keyword{ts}
