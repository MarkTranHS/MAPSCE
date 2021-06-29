#' index_tree
#'
#' `index_tree` indexes all branches in a given tree
#'
#' @param tree
#' @return tibble with indices
#'
#'
#' @keywords internal
#'

index_tree <- function(tree) {
  tree_nodes <- dplyr::tibble(node_id = unique(c(tree$parent, tree$child))) %>%
    dplyr::left_join(tree %>%
                       dplyr::rename(parent_id = parent), by = c("node_id" = "child")) %>%
    dplyr::mutate(left_index = NA, right_index = NA)
  root_id <- setdiff(unique(tree$parent), unique(tree$child))
  if (length(root_id) != 1) {
    stop("tree length is 1")
  }
  recursive_index <- function(tree_nodes, this_node_id) {
    index <- max(c(tree_nodes$left_index, tree_nodes$right_index, 0), na.rm = T) + 1
    tree_nodes <- tree_nodes %>%
      dplyr::mutate(left_index = ifelse(node_id == this_node_id, index, left_index))
    children_ids <- tree_nodes %>%
      dplyr::filter(parent_id == this_node_id) %>%
      dplyr::pull(node_id)
    for (this_child_id in children_ids) {
      tree_nodes <- recursive_index(tree_nodes, this_child_id)
    }
    index <- max(c(tree_nodes$left_index, tree_nodes$right_index, 0), na.rm = T) + 1
    tree_nodes <- tree_nodes %>%
      dplyr::mutate(right_index = ifelse(node_id == this_node_id, index, right_index))
    return(tree_nodes)
  }
  tree_nodes <- recursive_index(tree_nodes, root_id)
  return(tree_nodes)
}



