#' interpret_mapsce
#'
#' `interpret_mapsce` is a wrapper function for MAPSCE for interpretation of results for a particular mapping result
#'
#' @param data a tibble with summarised mapping results, the output of using mapsce function
#' @param tree a matrix listing all the branches in the tree, where the first column is the ancestral node and the second column is the descendant clone. By definition, the root node will only be present in the first column. The clone IDs must correspond to the cluster IDs in cluster_ccf and mutation_ccf.
#' @param min_diff a numeric vector, minimum difference between before and after to consider clonality. default at 0.4
#' @param consensus_threshold a numeric vector, the allowed threshold between the values to be considered "in agreement/consensus". default at 0.1
#' @param format a character vector, "list" or "tibble" output format. default is tibble
#' @param graph logical, graphical output
#'
#' @return returns interpretation of MAPSCE in the form of vector
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
#' interpret_mapsce(mapsce_result, example_tree)
#' # returns an interpretation vector for this example mapsce result
#'
#'
#' @export

## Interpretation function
interpret_mapsce <- function(data, tree, min_diff = 0.4, consensus_threshold = 0.1, format = "tibble", graph = F){
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

  data <- data %>% dplyr::mutate(good_result = dplyr::case_when(
    good_result == 1 & abs(before-after) >= min_diff ~ 1,
    good_result == 1 & abs(before-after) < min_diff ~ 0,
    T ~ 0
    ))

  if(all(data$good_result == 0)){
    data <- data %>% dplyr::mutate(good_result = ifelse(null == "yes", 2, 0))
  }

  consensus_mapping <- get_consensus(data, tree, consensus_threshold)
  number_of_results <-  data %>% dplyr::filter(good_result >= 1) %>% nrow()
  nclones <- nrow(data)
  consistency_filter <- F
  data <- data %>% dplyr::arrange(desc(good_result))

  # more than 1 good result
  if(number_of_results > 1){
    number_of_nonNAs <- sum(!is.na(consensus_mapping["consensus",]))
    number_of_NAs <- nclones - number_of_nonNAs
    if(number_of_nonNAs <= 1 & number_of_results - 1 < number_of_NAs){
    data[2:nrow(data),"good_result"] <- 0
    consistency_filter <- T
    consensus_mapping <- get_consensus(data, tree, consensus_threshold)
    }

  }
    good_row <- data %>% dplyr::filter(good_result >= 1)
    good_branches <- good_row %>% dplyr::pull(branch)
    null_branch <- data %>% dplyr::filter(null == "yes") %>% dplyr::pull(branch)
    before <- good_row %>% dplyr::pull(before)
    after <- good_row %>% dplyr::pull(after)

    cn_state <- consensus_mapping["consensus",] %>% unique
    number_of_states <- length(cn_state)
    cn_state <- sort(cn_state, decreasing = ifelse(all(before < after), F,T))
    cn_diff <- abs(cn_state[1] - cn_state[2])
    if(is.na(cn_diff)){
      cn_diff <- 0
      cn_state[2] <- cn_state[1]
      if(number_of_states>1){
        cn_state[2] <- NA
        cn_diff <- NA
        }
      }
    clonality <- dplyr::case_when(
      cn_diff >= min_diff ~ "subclonal",
      is.na(cn_diff) ~ "subclonal",
      TRUE ~ "null"
    )
    clonality <- unique(clonality)
    if(clonality == "null"){
      branch <- null_branch
      } else {
        branch <- good_branches
        }

    interpretation_result <- tibble::tibble(branch = list(branch),
                                            clonality = clonality,
                                            consistency_filter = consistency_filter,
                                            node_ids = list(colnames(consensus_mapping)),
                                            consensus_mapping = list(consensus_mapping[nrow(consensus_mapping),]))
    if(format == "list"){
      interpretation_result <- list(branch = branch,
                                    clonality = clonality,
                                    consistency_filter = consistency_filter,
                                    node_ids = list(colnames(consensus_mapping)),
                                    consensus_mapping = consensus_mapping)

    }
    return(interpretation_result)
}

# Number of good results: 2
# Number of clones: 4
# Number of non-NAs: 1
# (1<2)
# Number of non-NAs < number of clones – number of good results
# Number of good results -1
# Non-NAs > 1, NAs = number of good results - 1

