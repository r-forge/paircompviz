hasse <- function(e,                       # adjacency matrix (with edge sizes)
                  v=NULL,                  # names of the vertices
                  elab="",                 # names of the edges
                  ecol="black",            # edge label color
                  ebg="gray",              # edge line color
                  vcol="black",            # vertex label color
                  vbg="white",             # vertex background color
                  vsize=1,                 # vertex size
                  fvlab=".",               # fake vertex label
                  fvcol="black",           # fake vertex label color
                  fvbg="white",            # fake vertex background color
                  fvsize=1,                # fake vertex size
                  febg="black",            # fake edge line color
                  fesize=1,                # fake edge line size
                  main=paste("Hasse Diagram of", deparse(substitute(e))),
                  compress=FALSE) {        # compress edges by creating fake vertices

  # validate e
  stopifnot(is.matrix(e))
  stopifnot(nrow(e) == ncol(e))
  stopifnot(colnames(e) == rownames(e))
  stopifnot(!is.na(e))
  stopifnot(e >= 0)

  # validate v
  if (is.null(v)) {
    v <- colnames(e)
    if (is.null(v)) {
      v <- 1:ncol(e)
    }
  }
  stopifnot(is.vector(v))
  stopifnot(length(v) == ncol(e))
  stopifnot(!is.na(v))

  # validate elab
  if (is.null(elab)) {
    elab <- ""
  } 
  if (is.vector(elab) && (length(elab) == 1)) {
    tmp <- e
    tmp[tmp > 0] <- elab
    elab <- tmp
  } else {
    stopifnot(is.matrix(elab))
    stopifnot(nrow(elab) == ncol(elab))
    stopifnot(colnames(elab) == rownames(elab))
    stopifnot(nrow(elab) == nrow(e))
    stopifnot(colnames(elab) == colnames(e))
    stopifnot(!is.na(elab))
  }

  # validate ecol
  if (is.null(ecol)) {
    ecol <- "black"
  }
  if (is.vector(ecol) && (length(ecol) == 1)) {
    tmp <- e
    tmp[tmp > 0] <- ecol
    ecol <- tmp
  } else {
    stopifnot(is.matrix(ecol))
    stopifnot(nrow(ecol) == ncol(ecol))
    stopifnot(colnames(ecol) == rownames(ecol))
    stopifnot(nrow(ecol) == nrow(e))
    stopifnot(colnames(ecol) == colnames(e))
    stopifnot(!is.na(ecol))
  }

  # validate ebg
  if (is.null(ebg)) {
    ebg <- "gray"
  }
  if (is.vector(ebg) && (length(ebg) == 1)) {
    tmp <- e
    tmp[tmp > 0] <- ebg
    ebg <- tmp
  } else {
    stopifnot(is.matrix(ebg))
    stopifnot(nrow(ebg) == ncol(ebg))
    stopifnot(colnames(ebg) == rownames(ebg))
    stopifnot(nrow(ebg) == nrow(e))
    stopifnot(colnames(ebg) == colnames(e))
    stopifnot(!is.na(ebg))
  }

  # validate vcol
  if (is.null(vcol)) {
    vcol <- "black"
  }
  if (is.vector(vcol) && (length(vcol) == 1)) {
    vcol <- rep(vcol, length(v))
  } else {
    stopifnot(is.vector(vcol))
    stopifnot(length(vcol) == length(v))
    stopifnot(!is.na(vcol))
  }

  # validate vbg
  if (is.null(vbg)) {
    vbg <- "white"
  }
  if (is.vector(vbg) && (length(vbg) == 1)) {
    vbg <- rep(vbg, length(v))
  } else {
    stopifnot(is.vector(vbg))
    stopifnot(length(vbg) == length(v))
    stopifnot(!is.na(vbg))
  }

  # validate vsize
  if (is.null(vsize)) {
    vsize <- 1
  }
  if (is.vector(vsize) && (length(vsize) == 1)) {
    vsize <- rep(vsize, length(v))
  } else {
    stopifnot(is.vector(vsize))
    stopifnot(length(vsize) == length(v))
    stopifnot(!is.na(vsize))
    stopifnot(vsize > 0)
  }

  # validate fvlab
  if (is.null(fvlab)) {
    fvlab <- '.'
  }
  stopifnot(is.vector(fvlab))
  stopifnot(length(fvlab) == 1)

  # validate fvcol
  if (is.null(fvcol)) {
    fvcol <- 'black'
  }
  stopifnot(is.vector(fvcol))
  stopifnot(length(fvcol) == 1)

  # validate fvbg
  if (is.null(fvbg)) {
    fvbg <- 'white'
  }
  stopifnot(is.vector(fvbg))
  stopifnot(length(fvbg) == 1)

  # validate fvsize
  if (is.null(fvsize)) {
    fvsize <- 1
  }
  stopifnot(is.vector(fvsize))
  stopifnot(length(fvsize) == 1)
  stopifnot(fvsize > 0)

  # validate febg
  if (is.null(febg)) {
    febg <- 'black'
  }
  stopifnot(is.vector(febg))
  stopifnot(length(febg) == 1)

  # validate fesize
  if (is.null(fesize)) {
    fesize <- 1
  }
  stopifnot(is.vector(fesize))
  stopifnot(length(fesize) == 1)
  stopifnot(fesize > 0)

  stopifnot(is.logical(compress))


  # create graph
  graph <- as(e, "graphNEL")


  determineVertexName <- function(i) {
    if (i <= nrow(e)) {
      return(as.character(i))
    }
    return(paste("fake", i - nrow(e), sep=""))
  }


  # add fake edges if required
  fakeCount <- 0
  if (compress) {
    graphAdj <- adj(graph, nodes(graph))
    x <- 1
    while (x <= length(nodes(graph))) {
      y <- x+1
      while (y <= length(nodes(graph))) {
        intr <- intersect(graphAdj[[x]], graphAdj[[y]])
        if (length(intr) >= 2) {
          cc <- as.character(c(determineVertexName(x), determineVertexName(y)))
          z <- y+1
          while (z <= length(nodes(graph))) {
            if (length(intersect(graphAdj[[z]], intr)) == length(intr)) {
              cc <- as.character(c(cc, determineVertexName(z)))
            }
            z <- z+1
          }
          if (length(intr) * length(cc) >= length(intr) + length(cc)) {
            fakeCount <- fakeCount + 1
            fakeNode <- paste("fake", fakeCount, sep="")
            graph <- removeEdge(rep(cc, length(intr)),
                                rep(intr, times=rep(length(cc), length(intr))),
                                graph)
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


  # initialize node data (defaults for fake nodes and edges)
  nodeDataDefaults(graph, "shape") <- "rect"          # node shape
  nodeDataDefaults(graph, "size") <- fvsize           # node size
  nodeDataDefaults(graph, "fillcolor") <- fvbg        # node fill color
  nodeDataDefaults(graph, "textcolor") <- fvcol       # node font color
  edgeDataDefaults(graph, "lwd") <- fesize            # edge line width
  edgeDataDefaults(graph, "color") <- febg            # edge line color
  edgeDataDefaults(graph, "fontcolor") <- "black"     # edge label color
  edgeDataDefaults(graph, "label") <- ""              # edge label

  nodeData(graph, attr="shape") <- c(rep("rect", length(v)), rep("plaintext", fakeCount))
  nodeData(graph, attr="textcolor") <- c(vcol, rep(fvcol, fakeCount))
  nodeData(graph, attr="fillcolor") <- c(vbg, rep(fvbg, fakeCount))
  nodeData(graph, attr="size") <- c(vsize, rep(fvsize, fakeCount))

  # visualize pvalue by edge line width and color
  em <- t(edgeMatrix(graph))
  if (nrow(em) > 0) {
    for (x in 1:nrow(em)) {
      first <- em[x, 1]
      second <- em[x, 2]
      if (max(first, second) <= nrow(e)) {
        # NOT edge from/to fake node
        edgeData(graph, as.character(first), as.character(second), "lwd") <- e[first, second]
        edgeData(graph, as.character(first), as.character(second), "color") <- ebg[first, second]
        edgeData(graph, as.character(first), as.character(second), "fontcolor") <- ecol[first, second]
        edgeData(graph, as.character(first), as.character(second), "label") <- elab[first, second]
      }
    }
  }

  # default settings (will be applied mainly on fake nodes)
  attrs <- list(
    node = list(label=fvlab, shape="plaintext", width=1),
    edge = list(lwd=fesize, label="")
  )

  # pack everything
  nAttrs <- list()
  eAttrs <- list()
  nAttrs$shape <- nodeData(graph, attr="shape")
  nAttrs$fillcolor <- nodeData(graph, attr="fillcolor")
  nAttrs$fontcolor <- nodeData(graph, attr="textcolor")
  nAttrs$height <- nodeData(graph, attr="size")
  #nAttrs$width <- nodeData(graph, attr="size")

  # workaround to make the font smaller: if the labels are too short, add some spaces to them
  l <- v
  add <- pmax(0, ceiling((12 - nchar(l)) / 2))
  vectorizedRep <- Vectorize(function(x) { paste(rep(" ", x), collapse="") })
  l <- paste(vectorizedRep(add), l, vectorizedRep(add))
  names(l) <- c(as.character(1:length(v)))
  nAttrs$label <- l

  l <- edgeData(graph, attr="label")
  names(l) <- edgeNames(graph)
  eAttrs$label <- l

  l <- edgeData(graph, attr="fontcolor")
  names(l) <- edgeNames(graph)
  eAttrs$fontcolor <- l

  l <- edgeData(graph, attr="color")
  names(l) <- edgeNames(graph)
  eAttrs$color <- l

  l <- edgeData(graph, attr="lwd")
  names(l) <- edgeNames(graph)
  eAttrs$lwd <- l

  # plot the resulting graph
  if (main != "") {
    nrows <- length(which(strsplit(main, '')[[1]]=='\n')) + 1
    oldpar <- par(mar=c(0, 0, 3, 0), oma=c(0, 0, nrows, 0))
  }
  plot(graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs)
  if (main != "") {
    title(main=main, outer=TRUE)
    par(oldpar)
  }
}

