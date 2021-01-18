#' Plot a precision-recall curve for a list of protein-protein interactions
#' 
#' Assess the intrinsic quality of a CF-MS dataset by evaluating the 
#' precision of the interaction network for any number of interactions scored. 
#' The x-axis in the output plot represents the total number of interactions 
#' (and can optionally be subset to a limited range), while the y-axis shows
#' the precision for each point in the ranked list. 
#' 
#' Note that the x-axis
#' shows recall in terms of the total number of interactions in the CF-MS 
#' dataset, not the proportion of true-positive interactions recovered from the
#' 'gold standard' dataset, since many of these interactions may not have been
#' detected by CF-MS. 
#' 
#' @param pairs a data frame in which each row represents a protein pair, and
#'   which contains the following columns:
#' \enumerate{
#'   \item \code{score}: the score assigned to that protein pair, e.g., by a 
#'     machine-learning classifier, in which higher scores represent a greater 
#'     probability of a physical interaction
#'   \item \code{label}: a ground-truth annotation of whether that pair of
#'     proteins is known to physically interact; one of \code{1}, \code{0}, or 
#'     \code{NA} (not labelled) 
#' }
#' @param max_n optionally, plot the precision-recall curve only out to this
#'   maximum number of interactions
#' 
#' @return a \code{ggplot2} object plotting the precision-recall curve for the input 
#'   interactions
#' 
#' @importFrom tidyr drop_na
#' @importFrom AUC roc
#' @import ggplot2
#' @importFrom magrittr %>%
#' 
#' @export
plot_pr = function(pairs, max_n = NULL) {
  # arrange by score in descending order, if this was not already done
  pairs %<>% arrange(desc(score))
  
  # subset interactions
  pairs %<>%
    mutate(idx = row_number())
  if (!is.null(max_n)) {
    pairs %<>%
      filter(idx <= max_n)
  } 
  
  # calculate precision
  pairs %<>% mutate(precision = calculate_precision(label))
  
  # return a basic PR plot
  p = pairs %>%
    ggplot(aes(x = idx, y = precision)) +
    geom_path() +
    scale_x_continuous('# of interactions') +
    scale_y_continuous('Precision') +
    theme_bw()
  return(p)
}
