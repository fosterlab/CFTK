#' Evaluate the recovery of known protein complexes
#' 
#' Assess the intrinsic quality of a CF-MS dataset by evaluating its ability
#' to recover known protein complexes, using receiver operating characteristic
#' (ROC) analysis. In this analysis, the correlation coefficients between
#' every pair of proteins in the dataset are ranked, and compared to a binary
#' outcome variable reflecting whether the two proteins are in the same complex.
#' The area under the curve (AUC) is returned as a measure of the ability of the
#' CF-MS data to recover known complexes. This measure ranges from 0 to 1, 
#' with 1 representing perfect recovery, and 0.5 representing random recovery. 
#' 
#' @param pairs a matrix of dimensions (# of proteins) x (# of proteins), 
#'   scoring every possible protein pair, in which higher values reflect more 
#'   similar pairs, e.g. as returned by \link{score_pairs}
#' @param adj an adjacency matrix between all complex proteins, with 
#'   intra-complex pairs as \code{1}s and inter-complex pairs as \code{0}s,
#'   e.g. as returned by \link{to_adjacency_matrix}
#' 
#' @return the area under the receiver operating characteristic curve 
#' 
#' @importFrom reshape2 melt
#' @importFrom dplyr filter n_distinct
#' @importFrom tidyr drop_na
#' @importFrom AUC auc roc 
calculate_complex_auc = function(pairs, adj) {
  # convert pairs to a data frame
  pairs = reshape2::melt(pairs, varnames = c('protein1', 'protein2'),
                         value.name = 'cor') %>%
    drop_na() %>%
    filter(as.integer(protein1) < as.integer(protein2))
  
  # filter to proteins in the complexes
  pairs0 = pairs %>%
    filter(gene1 %in% rownames(adj), gene2 %in% rownames(adj))
  
  # calculate AUC
  x = pairs0$cor
  y = adj[as.matrix(pairs0[, 1:2])]
  if (nrow(pairs0) == 0) 
    stop("no proteins overlap with the complexes")
  if (n_distinct(y) < 2)
    stop("could not find two distinct labels")
  auroc = auc(roc(predictions = x,
                  labels = factor(y, levels = c('0', '1'))))
  return(auroc)
}
