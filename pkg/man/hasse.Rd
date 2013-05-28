\name{hasse}
\alias{hasse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualization of Hasse diagram specified by an adjacency matrix}
\description{
  Given an adjacency matrix, this function displays the corresponding Hasse diagram.
  This is a wrapper function for graph creation using the \pkg{Rgraphviz} package.
}
\usage{
hasse(e,
      v=NULL,
      elab="",
      ecol="black",
      ebg="gray",
      vcol="black",
      vbg="white",
      vsize=1,
      fvlab=".",
      fvcol="black",
      fvbg="white",
      fvsize=1,
      febg="black",
      fesize=1,
      main=paste("Hasse Diagram of", deparse(substitute(e))),
      compress=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{An adjacency matrix, with \eqn{e_{i,j}} indicating the edge size between vertices \eqn{i}
  and \eqn{j} (\eqn{e_{i,j} = 0} means no edge between \eqn{i} and \eqn{j}). The matrix must be
  rectangular with non-negative non-missing values.
}
  \item{v}{Vector of names of the vertices. If null, the vertex names will be obtained from column
  names of adjacency matrix \eqn{e}.
}
  \item{elab}{Labels of the edges. If it is a scalar value, all edges would have the same label.
  Otherwise, \eqn{elab} must be a rectangular matrix (similar to adjacency matrix \eqn{e}).
  A value on \eqn{i}-th row and \eqn{j}-th column is a label of the edge between vertex \eqn{i} and
  vertex \eqn{j}.
}
  \item{ecol}{Edge label color. If scalar, all edge labels have the same color. Otherwise,
  \eqn{ecol} must be in the form of adjacency matrix: a value on \eqn{i}-th row and \eqn{j}-th
  column is a color of the label of the edge between vertex \eqn{i} and vertex \eqn{j}.
}
  \item{ebg}{Edge line color. If scalar, all edges have the same color. Otherwise,
  \eqn{ebg} must be in the form of adjacency matrix: a value on \eqn{i}-th row and \eqn{j}-th
  column is a color of the edge between vertex \eqn{i} and vertex \eqn{j}.
}
  \item{vcol}{Vertex label color. If scalar, all vertices have the same label color. Otherwise,
  \eqn{vcol} must be a vector of the size corresponding to the number of vertices.
}
  \item{vbg}{Vertex background color. If scalar, all vertices have the same background color.
  Otherwise, \eqn{vcol} must be a vector of the size corresponding to the number of vertices.
}
  \item{vsize}{Vertex sizes. If scalar, all vertices have the same size in the image.
  Otherwise, \eqn{vsize} must be a vector of the size corresponding to the number of vertices.
}
  \item{fvlab}{Labels of "dot" vertices. Must be scalar.
}
  \item{fvcol}{"dot" vertex label color. Must be scalar.
}
  \item{fvbg}{"dot" vertex background color. Must be scalar.
}
  \item{fvsize}{"dot" vertex size. Must be scalar.
}
  \item{febg}{Color of edges introduced by edge compression. Must be scalar.
}
  \item{fesize}{Thickness of edges introduced by edge compression. Must be scalar and non-negative.
}
  \item{main}{Main title of the diagram.
}
  \item{compress}{
    \code{TRUE} if the edges should be compressed, i.e. if the maximum bi-cliques have to be found
    in the graph and replaced with a "dot" vertex. (See examples.)
}
}
\details{
  This function depicts a Hasse diagram specified with an adjacency matrix \eqn{e}.
  Hasse diagram is a visualization of partially ordered set, by drawing its transitive reduction as
  an oriented graph. Each vertex corresponds to an element of the set. There is an edge between
  vertex \eqn{i} and vertex \eqn{j} iff \eqn{i < j} and there is no \eqn{z} such that \eqn{i < z < j}.
  
  The function is also capable of edge compression via introducing the "dot" edges:
  Let \eqn{U}, \eqn{V} be two disjoint non-empty sets of edges,
  such that for each \eqn{u} from \eqn{U} and \eqn{v} from \eqn{V}, there exists an edge
  from \eqn{u} to \eqn{v}. (The number of such edges equals \eqn{|U| \cdot |V|}.) Starting from
  \eqn{|U| > 2} and \eqn{|V| > 2}, the Hasse diagram may become too complicated and hence confusing.
  Therefore a \code{compress} argument exists in this function that enables ``compression'' of
  the edges in such a way that a new ``dot'' node \eqn{w} is introduced and \eqn{|U| \cdot |V|}
  edges between sets \eqn{U} and \eqn{V} are replaced with \eqn{|U|+|V|} edges from set \eqn{U}
  to node \eqn{w} and from node \eqn{w} to set \eqn{V}.
}
\value{
  Nothing.
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
  \code{\link{paircomp}}
}
\examples{
  # linear order
  e <- matrix(c(0, 1, 1, 0, 0, 1, 0, 0, 0), nrow=3, byrow=TRUE)
  hasse(e)

  # prepare adjacency matrix
  m <- matrix(0, byrow=TRUE, nrow=5, ncol=5)
  m[3, 1] <- 1
  m[3, 2] <- 1
  m[4, 1] <- 9
  m[4, 2] <- 1
  m[5, 1] <- 1
  m[5, 2] <- 1
  m

  mc <- m
  mc[mc > 0] <- "red"
  ms <- m
  ms[ms > 0] <- "blue"

  # view m with default settings
  hasse(m, ebg="black")

  # view m WITHOUT edge compression and some fancy adjustments
  hasse(v=c("a", "b", "c", "d", "e"), 
             vcol=c(gray(0.5), gray(1), rep(gray(0), 3)), 
             vbg=gray(5:1/5), vsize=1:5, e=m, ecol=mc, ebg=ms, elab=m,
             compress=FALSE)

  # view m WITH edge compression and some fancy adjustments
  hasse(v=c("a", "b", "c", "d", "e"), 
             vcol=c(gray(0.5), gray(1), rep(gray(0), 3)), 
             vbg=gray(5:1/5), vsize=1:5, e=m, ecol=mc, ebg=ms, elab=m,
             compress=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
