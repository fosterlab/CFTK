#' Make labels for a list of protein pairs based on a gold standard
#' 
#' Create labels for a data frame containing protein pairs, by cross-referencing
#' those protein pairs with a 'gold standard' set of known interactions. 
#' The gold standard interactions can be provided in a variety of different
#' formats.
#' 
#' @param gold_standard a set of 'gold-standard' interactions. Can be provided
#'   in the form of an adjacency matrix, a data frame of protein pairs, or
#'   a list of protein complexes
#' @param dat a data frame with interacting proteins in the first two
#'   columns 
#' 
#' @return a vector of the same length as the number of rows in \code{dat}, 
#'   containing \code{NA}s for protein pairs not in the gold standard, and
#'   ones or zeroes otherwise
#' 
#' @importFrom tester is_square_matrix
#' 
#' @export
label_pairs = function(gold_standard, dat) {
  # convert gold standard to adjacency matrix, if it is not
  if (is.data.frame(gold_standard)) {
    ## pairwise data frame
    message("interpreting gold_standard as a data frame of protein pairs ...")
    adj = to_adjacency_matrix(gold_standard)
  } else if (is.list(gold_standard)) {
    message("interpreting gold_standard as a list of protein complexes ...")
    pairs = to_pairwise_df(gold_standard)
    adj = to_adjacency_matrix(pairs)
  } else if (is_square_matrix(gold_standard)) {
    message("interpreting gold_standard as an adjacency matrix ...")
    adj = gold_standard
  } else {
    stop("could not parse input gold_standard")
  }
  
  # label protein pairs
  proteins_1 = dat[[1]]
  proteins_2 = dat[[2]]
  lab_idxs = proteins_1 %in% rownames(adj) &
    proteins_2 %in% rownames(adj)
  if (sum(lab_idxs) == 0)
    stop("no proteins overlap between input and gold standard")
  idxing_mat = cbind(proteins_1[lab_idxs], proteins_2[lab_idxs])
  labels = rep(NA, nrow(dat))
  labels[lab_idxs] = adj[idxing_mat]
  
  return(labels)
}
