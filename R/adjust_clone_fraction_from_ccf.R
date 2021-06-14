#'adjust_clone_fraction_from_ccf
#'
#'`adjust_clone_fraction_from_ccf` returns adjusted clone CCF based on tree
#'
#'@param clone_fraction clone CCF
#'@param tree tree
#'@param node node
#'
#'@return adjusted clone CCF
#'
#' @keywords internal
#'

#Adjusting CCF to match the tree branches if there are less branches than trees
adjust_clone_fraction_from_ccf <- function(clone_fraction, tree, node) {
  for (descendent in tree[tree[, 1] == node, 2]) {
    clone_fraction[node, ] <- clone_fraction[node, ] - clone_fraction[descendent, ]
    clone_fraction <- adjust_clone_fraction_from_ccf(clone_fraction, tree, descendent)
  }
  return(clone_fraction)
}
