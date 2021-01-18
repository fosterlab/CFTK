#' Evaluate the recovery of proteins annotated to the same GO term
#' 
#' Assess the intrinsic quality of a CF-MS dataset by evaluating its ability
#' to recover groups of proteins annotated to the same GO term, using receiver
#' operating characteristic (ROC) analysis. In this analysis, the correlation 
#' coefficients between every pair of proteins in the dataset are ranked, and 
#' compared to a binary outcome variable reflecting whether the two proteins 
#' annotated to the same GO term. The analysis is then repeated for each GO term
#' in turn. The area under the curve (AUC) for each GO term is returned as a 
#' measure of the ability of the CF-MS data to recover proteins annotated to
#' this term. This measure ranges from 0 to 1, with 1 representing perfect 
#' recovery, and 0.5 representing random recovery. 
#' 
#' @param pairs a matrix of dimensions (# of proteins) x (# of proteins), 
#'   scoring every possible protein pair, in which higher values reflect more 
#'   similar pairs, e.g. as returned by \link{score_pairs}. Alternatively,
#'   a data frame of candidate protein-protein interactions, with proteins
#'   in the first two columns.
#' @param ann a list in which each entry corresponds to a GO term and contains
#'   all of the proteins annotated to that GO term, e.g. as returned by 
#'   \link{as_annotation_list}
#' @param score_column when \code{pairs} is a data frame, the column that 
#'   contains the score for each protein pair
#' @param verbose set to \code{FALSE} to disable messages from the function
#' 
#' @return a data frame with four columns:
#' \enumerate{
#'   \item \code{go_term}: the GO term in question, obtained from the names of 
#'     the input annotation list
#'   \item \code{n_proteins}: the number of proteins in the input annotation 
#'     list which are annotated to that GO term
#'   \item \code{n_chromatograms}: the number of proteins in the CF-MS dataset
#'     which are annotated to that GO term
#'   \item \code{auroc}: the AUC for that GO term
#' }
#'
#' 
#' @importFrom reshape2 melt
#' @importFrom dplyr filter mutate bind_rows
#' @importFrom tidyr drop_na
#' @importFrom AUC auc roc
#' @importFrom magrittr %<>%
calculate_go_auc = function(pairs, ann, score_column = 'cor', verbose = TRUE) {
  # convert pairs to a data frame
  if (!is.data.frame(pairs)) {
    pairs = reshape2::melt(pairs, varnames = c('protein1', 'protein2'),
                           value.name = 'cor') %>%
      drop_na() %>%
      filter(as.integer(protein1) < as.integer(protein2))
  } else {
    colnames(pairs)[c(1, 2)] = c('protein1', 'protein2')
  }
  
  # keep only proteins that have at least one GO annotation
  go_terms = unique(unlist(ann))
  proteins = with(pairs, unique(c(protein1, protein2)))
  keep = intersect(proteins, unlist(ann))
  pairs %<>% filter(protein1 %in% keep, protein2 %in% keep)
  
  # map over GO terms
  results = data.frame()
  for (term_idx in seq_along(ann)) {
    go_term = names(ann)[term_idx]
    if (verbose)
      message("[", term_idx, "/", length(ann), "] analyzing GO term: '",
              go_term, "' ...")
    targets = ann[[term_idx]]
    true_positives = intersect(targets, proteins)
    
    # create two vectors
    pairs0 = pairs %>%
      mutate(y = ifelse(protein1 %in% targets & protein2 %in% targets, 1, 0))
    x = pairs0$cor
    y = pairs0$y
    if (n_distinct(y) < 2) 
      next
    auroc = auc(roc(predictions = x,
                    labels = factor(y, levels = c('0', '1'))))
    
    # append results
    row = data.frame(go_term = go_term,
                     n_proteins = length(targets),
                     n_chromatograms = length(true_positives),
                     auroc = auroc)
    results %<>% bind_rows(row)
  }
  
  return(results)
}
