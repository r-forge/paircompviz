paircomp <- function(obj, 
                     grouping=NULL,
                     test=c("t", "prop", "wilcox"),
                     level=0.05,
                     main=NULL,
                     compress=TRUE,
                     visualize=c("position", "size", "pvalue"),
                     result=FALSE,
                     draw=TRUE,
                     ...) {

  stopifnot(is.logical(compress))
  stopifnot(is.logical(result))
  stopifnot(is.numeric(level))
  stopifnot(length(level) == 1)

  visualize <- match.arg(visualize, several.ok=TRUE)

  if (is.null(grouping)) {
    # obj must be glht from the multcomp package
    stopifnot(require("multcomp", quietly=TRUE))
    stopifnot(require("reshape", quietly=TRUE))
    stopifnot(is(obj, "glht"))
    stopifnot(length(obj$focus) == 1)

    testResult <- summary(obj)
    lev <- testResult$model$xlevels[[1]]
    pval <- as.data.frame(t(combn(lev, 2)))
    pval$pvals <- testResult$test$pvalues

    pval <- cast(pval, V2 ~ V1, value="pvals")
    pval <- as.matrix(pval)
    pval <- pval[lev[-1], lev[-length(lev)]]

    means <- tapply(testResult$model$model[, 1], testResult$model$model[, 2], mean)
    sizes <- tapply(testResult$model$model[, 1], testResult$model$model[, 2], length)

    if (is.null(main)) {
      objName <- deparse(substitute(obj))
      groupingName <- deparse(substitute(grouping))
      main <- paste("Pairwise comparisons on",
                    paste(colnames(testResult$model$model), collapse=" by "),
                    "\nusing",
                    sub('.', ' ', testResult$alternative, fixed=TRUE),
                    testResult$type,
                    "test")
    }

  } else {
    # perform pairwise tests
    objName <- deparse(substitute(obj))
    groupingName <- deparse(substitute(grouping))

    if (!is.factor(grouping)) {
      grouping <- factor(grouping)
    }

    stopifnot(is.vector(obj))
    test <- match.arg(test)

    nonnull <- !is.na(obj) & !is.na(grouping)
    obj <- obj[nonnull]
    grouping <- grouping[nonnull]

    if (test == "t") {
      means <- tapply(obj, grouping, mean)
      testResult <- pairwise.t.test(obj, grouping, ...)

    } else if (test == "wilcox") {
      means <- tapply(obj, grouping, median)
      testResult <- pairwise.wilcox.test(obj, grouping, ...)
      
    } else if (test == "prop") {
      positive <- tapply(obj, grouping, sum)
      count <- tapply(obj, grouping, length)
      means <- positive / count
      testResult <- pairwise.prop.test(positive, count, ...)
    }

    pval <- testResult$p.value
    lev <- levels(grouping)
    sizes <- tapply(obj, grouping, length)

    if (is.null(main)) {
      main <- paste("Pairwise comparisons on",
                    objName,
                    "by",
                    groupingName,
                    "\nusing",
                    testResult$method,
                    "and",
                    testResult$p.adjust.method,
                    "p-value adjustment")
    }
  }


  means <- as.vector(means)
  sizes <- as.vector(sizes)


  # create rectangular matrix of p-values
  pval <- rbind(NA, pval)
  pval <- cbind(pval, NA)
  rownames(pval) <- NULL
  colnames(pval) <- NULL
  pval <- ifelse(is.na(pval), t(pval), pval)
  pval <- ifelse(is.na(pval), 1, pval)

  # create adjacency matrix
  vertexCount <- length(means)
  x <- matrix(means, nrow=vertexCount, ncol=vertexCount)
  y <- matrix(means, nrow=vertexCount, ncol=vertexCount)

  strength <- pmax(pmin(floor(level / pval), 7), 0)
  adjw <- ifelse(t(x) < y, strength, 0)
  #adjw <- ifelse(t(x) < y,
                 #ifelse(pval <= 0.001, 7,
                 #ifelse(pval <= 0.005, 5,
                 #ifelse(pval <= 0.01, 3,
                 #ifelse(pval <= 0.05, 1, 0)))), 0)
  adjb <- transReduct(adjw)

  if (min(adjb) < 0) {
    stop("Broken transitivity -- cannot create Hasse diagram")
  }


  if ("position" %in% visualize && max(means) != min(means)) {
    bg=gray(0:100/100)
    fg=c(rep("white", floor(length(bg) / 2)), rep("black", ceiling(length(bg) / 2)))
    index <- 1 + round((length(bg) - 1) * (means - min(means)) / (max(means) - min(means)))
    vcol <- fg[index]
    vbg <- bg[index]
  } else {
    vcol <- "black"
    vbg <- "white"
  }


  if ("size" %in% visualize && max(sizes) > 0) {
    #vsize <- sqrt(sizes / max(sizes)) 
    vsize <- pmax(sizes / max(sizes), 0.2)
  } else {
    vsize <- 1
  }

  if ("pvalue" %in% visualize) {
    elab <- format(pval, digits=4)
    adj <- ifelse(adjb > 0, adjw, 0)
  } else {
    elab <- ""
    adj <- adjb
  }


  if (draw) {
    hasse(e=adj,
          v=lev,
          elab=elab,
          vsize=vsize,
          vcol=vcol,
          vbg=vbg,
          main=main,
          compress=compress)
  }
  if (result) {
    return(testResult)
  }
}
