\name{brokentrans}
\docType{data}
\alias{brokentrans}
\title{Artificial dataset that suffers with broken transitivity of the pairwise t-test comparisons}
\description{
  This is a dataset of artificial data created to demonstrate that there exists a data sample such
  that pairwise comparisons using t-test break transitivity of the results, i.e. that if treatment 1
  is significantly lower than 2 and treatment 2 is lower than 3, it is not always the case that also
  treatment 1 is significantly lower than 3.
}
\usage{data(brokentrans)}
\format{A data frame of two columns:
  \enumerate{
    \item x is the measured value,
    \item g is the treatment group.
  }
}
%\source{World Almanac and Book of Facts, 1975, page 406.}
%\references{
  %McNeil, D. R. (1977) \emph{Interactive Data Analysis}.
  %New York: Wiley.
%}
\examples{
  data(brokentrans)

  # For \alpha = 10^{-9}, we obtain significant difference
  # between 1-2, 2-3, but not 1-3.
  tapply(brokentrans$x, brokentrans$g, mean)
  pairwise.t.test(brokentrans$x, brokentrans$g, pool.sd=FALSE)
}
\keyword{datasets}
