\name{paircomp}
\alias{paircomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualization of multiple pairwise comparison test results}
\description{
  This function performs multiple pairwise comparison tests on given data and views the results in
  the form of Hasse diagram.
}
\usage{
paircomp(obj, 
         grouping=NULL,
         test=c("t", "prop", "wilcox"),
         level=0.05,
         main=NULL,
         compress=TRUE,
         visualize=c("position", "size", "pvalue"),
         result=FALSE,
         draw=TRUE,
         ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    either a vector or an object of class \code{\link[multcomp]{glht}} as returned from the
    \code{\link[multcomp]{glht}} function of the \code{multcomp} package.

    If \code{obj} is an object of class \code{\link[multcomp]{glht}}, then arguments \code{grouping}
    and \code{test} may be arbitrary, because they will be not used. Otherwise, if \code{test}
    equals \code{"t"} or \code{"wilcox"}, \code{obj} should be a numeric vector of responses, and
    if \code{test} equals \code{"prop"}, \code{obj} should be a vector of 0's and 1's.
}
  \item{grouping}{
    a grouping factor. If \code{obj} is a numeric vector, \code{grouping} must be a factor.
    If \code{obj} is an object of class \code{\link[multcomp]{glht}}, \code{grouping} should be
    \code{NULL}.
}
  \item{test}{
    a name of the test to use. If \code{obj} is an object of class \code{\link[multcomp]{glht}}, the value of
    \code{test} does not have any effect. Otherwise, the values determine the type of the pairwise
    comparison test procedure. Allowed values \code{"t"}, \code{"prop"} or \code{"wilcox"} imply
    internal call of \code{pairwise.t.test()}, \code{pairwise.prop.test()} or
    \code{pairwise.wilcox.test()}, respectively.
}
  \item{level}{
    the maximum p-value that will be considered as significant; i.e. pairwise test results with
    p-value lower than the specified level will be represented with an edge in the resulting Hasse
    diagram.
}
  \item{main}{
    main title of the diagram.
}
  \item{compress}{
    \code{TRUE} if the edges should be compressed, i.e. if the maximum bi-cliques have to be found
    in the graph and replaced with a "dot" vertex. (See examples.)
}
  \item{visualize}{
    vector of additional information to be included in the diagram: \code{"position"}
    enables vertex background color to be derived from the treatment's proportion (\code{"prop"}
    test) or mean value (otherwise); \code{"size"} enables vertex size to correspond to the
    treatment's sample size; \code{"pvalue"} sets the edges' line thickness accordingly to p-value
    (lower p-value corresponds to thicker line).
}
  \item{result}{
    whether to return test results as a return value.
}
  \item{draw}{
    whether to render the diagram.
}
  \item{...}{
    other arguments that will be passed to the underlying function that performs pairwise
    comparisons (e.g. \code{pairwise.t.test}, \code{pairwise.prop.test} or
    \code{pairwise.wilcox.test}.
}
}
\details{
  All treatments in a set are compared in pairs using a selected statistical test. If the results
  form a partially ordered set, they can be viewed in a Hasse diagram.
  
  Hasse diagram is a graph with each treatment being represented as a vertex. An edge is drawn
  downwards from vertex \eqn{y} to vertex \eqn{x} if and only if treatment \eqn{x} is significantly
  lower than treatment \eqn{y}, and there is no such treatment \eqn{z} that \eqn{x} was lower than
  \eqn{z} and \eqn{z} lower than \eqn{y}. Each edge is connected to exactly two vertices: its two
  endpoints. If there does not exist a path between some two treatments, it means that these two
  treatments are incomparable (i.e. the difference among them is not statistically significant).

  The function accepts two types of inputs: either an instance of class \code{\link[multcomp]{glht}}
  or a vector \code{obj} of measured values and a factor \code{grouping} of treatments. 
  
  The \code{\link[multcomp]{glht}} object may be obtained from function \code{\link[multcomp]{glht}} 
  of the \pkg{multcomp} package and set as the \code{obj} argument. Argument \code{grouping} must be
  \code{NULL}, in that case.

  If \code{obj} is a numeric vector of measured values, \code{grouping} must not be \code{NULL}
  and also a type of statistical test must be selected by setting \code{test} argument.

  Edge compression (introducing "dot" edges):

  Sometimes, pairwise comparison tests may yield in such bipartite setting that each pair of nodes
  of some two node subsets would be inter-connected with an edge (without any edge between nodes in
  the same subset). More specifically, let \eqn{U}, \eqn{V} be two disjoint non-empty sets of edges,
  such that for each \eqn{u} from \eqn{U} and \eqn{v} from \eqn{V}, there exists an edge
  from \eqn{u} to \eqn{v}. (The number of such edges equals \eqn{|U| \cdot |V|}.) Starting from
  \eqn{|U| > 2} and \eqn{|V| > 2}, the Hasse diagram may become too complicated and hence confusing.
  Therefore a \code{compress} argument exists in this function that enables ``compression'' of
  the edges in such a way that a new ``dot'' node \eqn{w} is introduced and \eqn{|U| \cdot |V|}
  edges between sets \eqn{U} and \eqn{V} are replaced with \eqn{|U|+|V|} edges from set \eqn{U}
  to node \eqn{w} and from node \eqn{w} to set \eqn{V}.
}
\value{
  If argument \code{result} is \code{TRUE}, the function returns everything that is returned by the
  underlying test function (\code{\link{pairwise.t.test}}, \code{\link{pairwise.prop.test}} or
  \code{\link{pairwise.wilcox.test}} accordingly to the \code{test} argument), or a copy of the
  \code{obj} argument, if \code{obj} is an instance of class \code{\link[multcomp]{glht}}.
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
  \code{\link{pairwise.t.test}},
  \code{\link{pairwise.prop.test}},
  \code{\link{pairwise.wilcox.test}},
  \code{\link[multcomp]{glht}}
  \code{\link{hasse}}
}
\examples{
  # Example of test="prop":
  o <- c(rep(1, 10), rep(0, 10), rep(c(0,1), 5))
  g <- c(rep(1,10), rep(2, 10), rep(3, 10))
  paircomp(o, g, test="prop")

  # Example of test="t" and test="wilcox":
  paircomp(InsectSprays$count, InsectSprays$spray, test="t")
  paircomp(InsectSprays$count, InsectSprays$spray, test="wilcox")

  # Example of t-test with non-pooled SD and Bonferroni adjustment
  # for multiple comparisons:
  paircomp(InsectSprays$count, InsectSprays$spray, test="t",
           pool.sd=FALSE, p.adjust.method="bonferroni")

  # Compare diagrams with and without compressed edges:
  paircomp(InsectSprays$count, InsectSprays$spray, test="t",
           compress=FALSE)
  paircomp(InsectSprays$count, InsectSprays$spray, test="t",
           compress=TRUE)

  # perform Tukey test:
  library(rpart) # for car90 dataset
  library(multcomp) # for glht() function
  aovR <- aov(Price ~ Type, data = car90)
  glhtR <- glht(aovR, linfct = mcp(Type = "Tukey"))
  paircomp(glhtR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ htest }
