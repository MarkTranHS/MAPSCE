#' get_consensus
#'
#' `get_consensus` finds a consensus copy number results for a particular mapping result
#'
#' @param data summarised mapping results from mapsce function
#' @param tree tree
#' @param consensus printing of the consensus
#' @param consensus.only limiting output to just the consensus
#' @return returns a matrix with consensus copy number states for every clone in a tree, labeling no consensus as NAs
#'
#' @examples
#' data(example_data)
#'
#' mapsce_result <- mapsce(
#' example_cn,
#' example_ccf,
#' example_mutational_ccf,
#' example_tree)
#'
#' get_consensus(mapsce_result, example_tree)
#' # returns a consensus matrix for this example mapsce result
#'
#'
#' @export

## Consensus function
get_consensus <- function(data, tree, consensus=TRUE, consensus.only=F){
  null <- data %>%
    dplyr::filter(null == "yes") %>%
    dplyr::pull(branch)
  all_branches <- data %>%
    dplyr::filter(evid >= 1) %>%
    dplyr::pull(branch) %>% unique

  tree_nodes <- NULL


  tree <- tibble::tibble(parent = tree[,1], child = tree[,2])

  tree_nodes <- tibble::tibble(node_id = c(1:length(unique(union(tree$parent, tree$child)))),
                                branch = tree$parent %>% union(tree$child) %>% unique,
                                left_index = NA,
                                right_index = NA)


  # indexing_branches <- function(this_node, this_index) {
  #   if (missing(this_node)) {
  #     this_node <- 1
  #   } else {
  #     this_node <- tree_nodes %>%
  #       dplyr::filter(branch==this_node) %>%
  #       dplyr::pull(node_id)
  #   }
  #   if (missing(this_index)) {
  #     this_index <- 1
  #   }
  #
  #
  #   tree_nodes <<- tree_nodes %>%
  #     dplyr::mutate(left_index = ifelse(node_id == this_node, this_index, left_index))
  #   this_index <- this_index +1
  #   node <- tree_nodes %>%
  #     dplyr::filter(node_id == this_node) %>%
  #     dplyr::pull(branch)
  #   children <- tree %>%
  #     dplyr::filter(parent == node) %>%
  #     dplyr::pull(child)
  #   for (this_child_node_id in children) {
  #     this_index <- indexing_branches(this_child_node_id, this_index)
  #   }
  #   tree_nodes <<- tree_nodes %>%
  #     dplyr::mutate(right_index = ifelse(node_id == this_node, this_index, right_index))
  #   this_index <- this_index + 1
  #   return(invisible(this_index))
  # }
  #indexing_branches()
  index_tree(tree)
  tree_nodes <- tree_nodes %>%
    dplyr::select(-node_id)

  c_tree <- matrix(nrow = length(all_branches), ncol = nrow(tree_nodes))
  colnames(c_tree) <- tree_nodes[,1] %>% dplyr::pull()
  rownames(c_tree) <- all_branches
  for(this_branch in all_branches){
    cn_children <- data %>% dplyr::filter(branch == this_branch) %>% dplyr::pull(after)
    cn_parents <- data %>% dplyr::filter(branch == this_branch) %>% dplyr::pull(before)
    left <- tree_nodes %>% dplyr::filter(branch == this_branch) %>% dplyr::pull(left_index)
    right <- tree_nodes %>% dplyr::filter(branch == this_branch) %>% dplyr::pull(right_index)
    children <- c(tree_nodes %>%
                    dplyr::filter(left_index > left & right_index < right) %>%
                    dplyr::pull(branch),this_branch)
    parents <- setdiff(tree_nodes %>% dplyr::pull(branch), children)
    for(child in children){
      c_tree[as.character(this_branch),child] <- round(cn_children,digits=2)
    }
    for(parent in parents){
      c_tree[as.character(this_branch),parent] <- round(cn_parents,digits=2)
    }
  }

  if(consensus==TRUE){
    consensus <- apply(c_tree, 2, function(x) {ifelse(all(abs(x - mean(x))<=0.1), mean(x), NA) } )
    c_tree <- rbind(c_tree,consensus)

    if(consensus.only==TRUE){
      c_tree <- consensus
    }
  }

  return(c_tree)
}

