#' Convert an annotations list to a pairwise data frame
#' 
#' Convert a list of gene sets (for instance, where each entry corresponds to a
#' complex and contains all the protein subunits of that complex) into a data
#' frame enumerating all of the possible pairs within each list item. 
#'
#' @param ann a list of gene sets, e.g. as returned by \link{as_annotation_list}
#' @return a data frame with two columns, \code{'protein_A'} and 
#'   \code{'protein_B'}, containing all unique pairs of protiens found within 
#'   the same gene set
#' 
#' @importFrom tidyr crossing
#' @importFrom dplyr filter distinct
#' @importFrom purrr map_dfr
#' @importFrom magrittr %>%
#' 
#' @export
to_pairwise_df = function(ann) {
  ann %>%
    map_dfr(~ tidyr::crossing(protein_A = ., protein_B = .), 
            .id = 'complex') %>%
    filter(protein_A < protein_B) %>%
    distinct(protein_A, protein_B)
}
