\name{hilbert.filter}
\alias{hilbert.filter}
\title{Select a Hilbert Wavelet Pair}
\description{
  Converts name of Hilbert wavelet pair to filter coefficients.
}
\usage{
hilbert.filter(name)
}
\arguments{
  \item{name}{Character string of Hilbert wavelet pair, see acceptable
    names below (e.g., \code{"k3l3"}).}
}
\value{
  List containing the following items:
  \item{L}{length of the wavelet filter}
  \item{h0,g0}{low-pass filter coefficients}
  \item{h1,g1}{high-pass filter coefficients}
}
\details{
  Simple \code{switch} statement selects the appropriate HWP.  There are
  two parameters that define a Hilbert wavelet pair using the notation
  of Selesnick (200X), $K$ and $L$.  Currently, the only implemented
  combinations $(K,L)$ are (3,3), (3,5), (4,2) and (4,4).
}
\references{
  Selesnick, I.W. (2001).  Hilbert transform pairs of wavelet bases. 
  \emph{IEEE Signal Processing Letters\/}~\bold{8}(6), 170--173.

  Selesnick, I.W. (2002).  The design of approximate {H}ilbert transform
  pairs of wavelet bases.  \emph{IEEE Transactions on Signal
    Processing\/}~\bold{50}(5), 1144--1152.
}
\seealso{
  \code{\link{wave.filter}}
}
%\examples{}
\author{B. Whitcher}
\keyword{ts}
