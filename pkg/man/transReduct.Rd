\name{transReduct}
\alias{transReduct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Remove transitive edges from an adjacency matrix}
\description{
  This function removes transitive edges from an adjacency matrix.
}
\usage{
transReduct(e)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{
    an adjacency matrix, i.e. a rectangular matrix with value \eqn{e_{i,j}} above zero indicating an
    edge between vertices \eqn{i} and \eqn{j} of the corresponding graph.
}
}
\details{
  This function takes an adjacency matrix as the argument \code{e}. Both rows and columns correspond
  to graph vertices, with value \eqn{e_{i,j}} above zero indicating an edge between vertices \eqn{i}
  and \eqn{j}. The function removes all transitive edges, i.e. sets to zero corresponding elements
  of matrix \eqn{e}. The transitive edge is such an edge between vertices \eqn{i} and \eqn{j} 
  that after removing it from the graph, there still exists a path from \eqn{i} to \eqn{j}.
}
\value{
  An adjacency matrix \eqn{e} with transitive edges being removed.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Michal Burda
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{paircomp}},
  \code{\link{hasse}}
}
\examples{
  e <- matrix(c(0, 1, 1, 0, 0, 1, 0, 0, 0), nrow=3, byrow=TRUE)
  transReduct(e)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
