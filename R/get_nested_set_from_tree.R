#' get_nested_set_from_tree
#'
#' `get_nested_set_from_tree` returns a nested set from a tree
#'
#' @param tree tree
#' @return a matrix
#' @keywords internal
#'


#Making a nested set matrix - 1/0
get_nested_set_from_tree <- function(tree) {
  clones <- unique(c(tree[, 1], tree[, 2]))
  number_of_clones <- length(clones)

  nested_set <- matrix(0, number_of_clones, number_of_clones)
  rownames(nested_set) <- clones
  colnames(nested_set) <- clones
  nested_set[tree] <- 1
  n <- 0
  while (n < sum(nested_set)) {
    n <- sum(nested_set)
    ancestors <- list()
    for (i in 1:number_of_clones) {
      ancestors[[i]] <- which(nested_set[, i] == 1)
    }
    for (i in 1:number_of_clones) {
      for (j in ancestors[[i]]) {
        nested_set[ancestors[[j]], i] <- 1
      }
    }
  }
  return(nested_set)
}
