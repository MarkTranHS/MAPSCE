#' reduce_pyclone_ccf
#'
#' `reduce_pyclone_ccf` returns a reduced cluster ccf
#'
#'@param ccf cluster ccf
#'@param tree tree
#'
#' @return reduced cluster ccf
#' @keywords internal
#'

#Reducing pyclone trees
reduce_pyclone_ccf <- function(ccf, tree) {
  nodes <- unique(c(tree[, 1], tree[, 2]))
  pyclone_ccf <- ccf[nodes, , drop = FALSE]
  return(pyclone_ccf)
}
