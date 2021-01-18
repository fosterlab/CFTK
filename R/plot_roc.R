#' Plot a ROC curve for a list of protein-protein interactions
#' 
#' Assess the intrinsic quality of a CF-MS dataset by evaluating its ability
#' to recover known protein complexes, using receiver operating characteristic
#' (ROC) analysis. This function visualizes the trade-off between the proportion 
#' of true-positive interactions recovered at any given rate of false-positive 
#' interactions. In the ideal analysis, this trade-off is skewed towards the top
#' left of the plot area, recovering all true positives without any false 
#' positives. The dashed line shows the performance of a theoretical random 
#' classifier. 
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
#' 
#' @return a \code{ggplot2} object plotting the ROC curve for the input 
#'   interactions
#' 
#' @importFrom tidyr drop_na
#' @importFrom AUC roc
#' @import ggplot2
#' @importFrom magrittr %>%
#' 
#' @export
plot_roc = function(pairs) {
  # remove unlabelled pairs
  pairs0 = drop_na(pairs, label)
  
  # perform ROC analysis
  x = pairs0$score
  y = factor(pairs0$label, levels = c('0', '1'))
  roc = roc(predictions = x, labels = y)
  
  # extract data frame
  df = data.frame(fpr = roc$fpr, tpr = roc$tpr)
  
  # return a basic ROC plot
  p = df %>%
    ggplot(aes(x = fpr, y = tpr)) +
    geom_abline(aes(slope = 1, intercept = 0), size = 0.4,
                linetype = 'dotted') +
    geom_path() +
    scale_x_continuous('False positive rate') +
    scale_y_continuous('True positive rate') +
    theme_bw() +
    theme(aspect.ratio = 1)
  return(p)
}
