\name{my.acf}
\alias{my.acf}
\title{Auotocovariance Function via the Discrete Fourier Transform}
\description{
  Computes the autocovariance function for a time series.
}
\usage{my.acf(x)
}
\arguments{
  \item{x}{time series}
}
\value{
  The autocovariance function, for all nonnegative lags, is output.
}
\details{
  The series is zero padded to twice its length before the discrete
  Fourier transform is applied.  Only the values corresponding to
  nonnegative lags are provided.
}
%\references{}
%\seealso{}
\examples{
data(ibm)
ibm.returns <- diff(log(ibm))
plot(1:length(ibm.returns) - 1, my.acf(ibm.returns), type="h",
     xlab="lag", ylab="ACVS", main="Autocovariance Sequence for IBM Returns")
}
\author{B. Whitcher}
\keyword{ts}
