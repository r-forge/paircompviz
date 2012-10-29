plot.pairwise.htest <- function(testResult,
                                data,
                                group,
                                compressEdges=TRUE,
                                fillColor=gray(0:100/100),
                                textColor=c(rep("white", floor(length(fillColor) / 2)), rep("black", ceiling(length(fillColor) / 2))),
                                fakeColorIndex=length(fillColor)) {
  stopifnot(inherits(testResult, "pairwise.htest"))
  stopifnot(is.vector(data))
  stopifnot(is.factor(group))
  stopifnot(length(fillColor) == length(textColor))
  stopifnot(length(fakeColorIndex) == 1)
  stopifnot(fakeColorIndex > 0)
  stopifnot(fakeColorIndex <= length(fillColor))

  vertexCount <- nlevels(group)
  means <- tapply(data, group, mean, na.rm = TRUE)

  # create rectangular matrix of p-values
  pval <- (testResult$p.value)
  pval <- rbind(NA, pval)
  pval <- cbind(pval, NA)
  rownames(pval) <- NULL
  colnames(pval) <- NULL
  pval <- ifelse(is.na(pval), t(pval), pval)
  pval <- ifelse(is.na(pval), 1, pval)

  # create adjacency matrix
  x <- matrix(means, nrow=vertexCount, ncol=vertexCount)
  y <- matrix(means, nrow=vertexCount, ncol=vertexCount)
  adjacencyWeighted <- ifelse(t(x) < y,
    ifelse(pval <= 0.001, 4,
    ifelse(pval <= 0.005, 3,
    ifelse(pval <= 0.01, 2,
    ifelse(pval <= 0.05, 1, 0)))), 0)
  adjacencyBinary <- ifelse(adjacencyWeighted > 0, 1, 0)

  # remove transitive edges
  transitivity <- ifelse(adjacencyBinary %*% adjacencyBinary > 0, 1, 0)

  adjacencyDirect <- adjacencyBinary - transitivity
  if (min(adjacencyDirect) < 0) {
    stop("Broken transitivity")
  }

  # create graph
  graph <- as(adjacencyDirect, "graphNEL")

  # add fake edges if required
  if (compressEdges) {
    fakeCount <- 0
    graphAdj <- adj(graph, nodes(graph))
    x <- 1
    while (x <= length(nodes(graph))) {
      y <- x+1
      while (y <= length(nodes(graph))) {
        intr <- intersect(graphAdj[[x]], graphAdj[[y]])
        if (length(intr) >= 2) {
          cc <- as.character(c(x, y))
          z <- y+1
          while (z <= length(nodes(graph))) {
            if (length(intersect(graphAdj[[z]], intr)) == length(intr)) {
              cc <- as.character(c(cc, z))
            }
            z <- z+1
          }
          if (length(intr) * length(cc) >= length(intr) + length(cc)) {
            fakeCount <- fakeCount + 1
            fakeNode <- paste("fake", fakeCount, sep="")
            graph <- removeEdge(rep(cc, length(intr)), rep(intr, times=rep(length(cc), length(intr))), graph)
            graph <- addNode(fakeNode, graph, edges=list(intr))
            graph <- addEdge(cc, fakeNode, graph, 1)
            graphAdj <- adj(graph, nodes(graph))
          }
        }
        y <- y+1
      }
      x <- x+1
    }
  }

  # initialize node data
  nodeDataDefaults(graph, "shape") <- "rect"       # node shape
  nodeDataDefaults(graph, "size") <- 1             # node size
  nodeDataDefaults(graph, "fillcolor") <- "white"  # node fill color
  nodeDataDefaults(graph, "textcolor") <- "black"  # node font color
  edgeDataDefaults(graph, "lwd") <- 1              # edge line width
  edgeDataDefaults(graph, "color") <- "black"      # edge line color

  # visualize means by node colors
  color <- 1 + round(length(fillColor) * (means - min(means)) / (max(means) - min(means)))
  color <- c(color, rep(fakeColorIndex, fakeCount))
  nodeData(graph, attr="fillcolor") <- fillColor[color]
  nodeData(graph, attr="textcolor") <- textColor[color]

  # change the shape
  nodeData(graph, attr="shape") <- c(rep("rect", vertexCount), rep("plaintext", fakeCount))

  # visualize pvalue by edge line width and color
  e <- t(edgeMatrix(graph))
  for (x in 1:nrow(e)) {
    first <- e[x, 1]
    second <- e[x, 2]
    if (max(first, second) < nrow(adjacencyWeighted)) {
      # NOT edge from/to fake node
      edgeData(graph, as.character(first), as.character(second), "lwd") <- (-1) + 1.5 * adjacencyWeighted[first, second]
      edgeData(graph, as.character(first), as.character(second), "color") <- gray(4:0/8)[(-1) + adjacencyWeighted[first, second]]
    }
  }

  # visualize variance by node size
  size <- sqrt(tapply(data, group, var, na.rm = TRUE))
  size <- 1 * size / max(size)
  size <- c(size, rep(min(size), fakeCount))
  nodeData(graph, attr="size") <- size

  # default settings (will be applied mainly on fake nodes)
  attrs <- list(
    node = list(label = ".", shape = "plaintext"),
    edge = list(lwd = 1)
  )

  # pack everything
  nAttrs <- list()
  eAttrs <- list()
  nAttrs$shape <- nodeData(graph, attr="shape")
  nAttrs$fillcolor <- nodeData(graph, attr="fillcolor")
  nAttrs$fontcolor <- nodeData(graph, attr="textcolor")
  nAttrs$height <- nodeData(graph, attr="size")
  nAttrs$width <- nodeData(graph, attr="size")
  eAttrs$lwd <- edgeData(graph, attr="lwd")

  l <- c(names(testResult$p.value[1,])[1], names(testResult$p.value[,1]))
  l <- c(levels(group))
  names(l) <- c(as.character(1:vertexCount))
  nAttrs$label <- l

  # for some strange reason I cannot do simply:
  #eAttrs$color <- edgeData(graph, attr="color")
  z <- edgeData(graph, attr="color")
  zz <- as.vector(z, mode="character")
  names(zz) <- sub("|", "~", names(z), fixed=TRUE)
  eAttrs$color <- zz

  # plot the resulting graph
  plot(graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs)
}

