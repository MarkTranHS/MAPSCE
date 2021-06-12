#' get_clone_fraction_from_ccf
#'
#' `get_clone_fraction_from_ccf` returns clone CCF from cluster CCF
#'
#' @param ccf ccf
#' @param tree tree
#'
#' @return a matrix
#' @keywords internal
#' @export
#'

#Setting trunk from tree to feed into recursive function
get_clone_fraction_from_ccf <- function(ccf, tree) {
  trunk <- setdiff(tree[, 1], tree[, 2])
  stopifnot(length(trunk) == 1)
  clone_fraction <- adjust_clone_fraction_from_ccf(ccf, tree, trunk)
  return(clone_fraction)
}
