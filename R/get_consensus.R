#' get_consensus
#'
#' `get_consensus` finds a consensus copy number results for a particular mapping result
#'
#' @param data a tibble with summarised mapping results, the output of using mapsce function
#' @param tree a matrix listing all the branches in the tree, where the first column is the ancestral node and the second column is the descendant clone. By definition, the root node will only be present in the first column. The clone IDs must correspond to the cluster IDs in cluster_ccf and mutation_ccf.
#' @param sens sensitivity of the comparison between the mapped values of different branches
#' @param consensus printing of the consensus
#' @param consensus.only limiting output to just the consensus, will output a vector instead of a data frame, only works with consensus = TRUE
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
get_consensus <- function(data, tree, sens = 0.1, consensus=TRUE, consensus.only=F){
  if(length(tree) == 0) {
    stop("missing tree")
  }
  if(is.data.frame(data) == F) {
    stop("needs data frame/tibble as input data")
  }
  if(!all(c("null", "branch", "good_result", "before", "after") %in% colnames(data))) {
    stop("function needs MAPSCE output with specific column names")
  }
  if(nrow(data) != length(unique(as.vector(tree)))) {
    stop("mismatch between data input and number of tree nodes")
  }
  null <- data %>%
    dplyr::filter(null == "yes") %>%
    dplyr::pull(branch)
  all_branches <- data %>%
    dplyr::filter(good_result >= 1) %>%
    dplyr::pull(branch) %>% unique

  tree_nodes <- NULL


  tree <- tibble::tibble(parent = tree[,1], child = tree[,2])

  tree_nodes <- tibble::tibble(node_id = c(1:length(unique(union(tree$parent, tree$child)))),
                                branch = tree$parent %>% union(tree$child) %>% unique,
                                left_index = NA,
                                right_index = NA)


  tree_nodes <- index_tree(tree) %>%
    dplyr::rename(branch = node_id)

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
    consensus <- apply(c_tree, 2, function(x) {ifelse(all(abs(x - mean(x))<= sens), mean(x), NA) } )
    c_tree <- rbind(c_tree,consensus)

    if(consensus.only==TRUE){
      c_tree <- consensus
    }
  }

  return(c_tree)
}

