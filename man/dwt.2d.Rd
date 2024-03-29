% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_D.R
\name{dwt.2d}
\alias{dwt.2d}
\alias{idwt.2d}
\title{Two-Dimensional Discrete Wavelet Transform}
\usage{
dwt.2d(x, wf, J = 4, boundary = "periodic")

idwt.2d(y)
}
\arguments{
\item{x}{input matrix (image)}

\item{wf}{name of the wavelet filter to use in the decomposition}

\item{J}{depth of the decomposition, must be a number less than or equal to
log(minM,N,2)}

\item{boundary}{only \code{"periodic"} is currently implemented}

\item{y}{an object of class \code{dwt.2d}}
}
\value{
List structure containing the \eqn{3J+1} sub-matrices from the
decomposition.
}
\description{
Performs a separable two-dimensional discrete wavelet transform (DWT) on a
matrix of dyadic dimensions.
}
\details{
See references.
}
\examples{

## Xbox image
data(xbox)
xbox.dwt <- dwt.2d(xbox, "haar", 3)
par(mfrow=c(1,1), pty="s")
plot.dwt.2d(xbox.dwt)
par(mfrow=c(2,2), pty="s")
image(1:dim(xbox)[1], 1:dim(xbox)[2], xbox, xlab="", ylab="",
      main="Original Image")
image(1:dim(xbox)[1], 1:dim(xbox)[2], idwt.2d(xbox.dwt), xlab="", ylab="",
      main="Wavelet Reconstruction")
image(1:dim(xbox)[1], 1:dim(xbox)[2], xbox - idwt.2d(xbox.dwt),
      xlab="", ylab="", main="Difference")

## Daubechies image
data(dau)
par(mfrow=c(1,1), pty="s")
image(dau, col=rainbow(128))
sum(dau^2)
dau.dwt <- dwt.2d(dau, "d4", 3)
plot.dwt.2d(dau.dwt)
sum(plot.dwt.2d(dau.dwt, plot=FALSE)^2)

}
\references{
Mallat, S. (1998) \emph{A Wavelet Tour of Signal Processing},
Academic Press.

Vetterli, M. and J. Kovacevic (1995) \emph{Wavelets and Subband Coding},
Prentice Hall.
}
\seealso{
\code{\link{modwt.2d}}.
}
\author{
B. Whitcher
}
\keyword{ts}
