\name{paircompviz-package}
\alias{paircompviz-package}
\alias{paircompviz}
\docType{package}
\title{Multiple comparison test visualization}
\description{
  This package provides visualization of the results from the multiple (i.e. pairwise)
  comparison tests such as pairwise.t.test, pairwise.prop.test or pairwise.wilcox.test. The groups
  being compared are visualized as nodes in Hasse diagram. Such approach enables very clear and
  vivid depiction of which group is significantly greater than which others, especially if comparing
  a large number of groups.
}
%\details{
%~~ An overview of how to use the package, including the most important functions ~~
%}
\author{
Michal Burda <michal.burda@osu.cz>
}
%\references{
%~~ Literature or other references for background information ~~
%}
\keyword{ hplot }
\keyword{ htest }
%\seealso{
  %\code{\link[multcomp:multcomp-package]{multcomp}} ~~
  %\code{\link[Rgraphviz:Rgraphviz-package]{Rgraphviz}} ~~
%}
\examples{
  paircomp(InsectSprays$count, InsectSprays$spray, test="t")
}
