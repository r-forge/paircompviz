transReduct <- function(e) {
  stopifnot(is.matrix(e))
  stopifnot(nrow(e) == ncol(e))
  stopifnot(!is.na(e))
  stopifnot(e >= 0)

  e <- ifelse(e > 0, 1, 0)
  transitivity <- ifelse(e %*% e > 0, 1, 0)
  return(e - transitivity)
}

