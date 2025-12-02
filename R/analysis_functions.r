#' Compute sensitivity, specificity and false positive rate across thresholds
#'
#' This function computes classification performance metrics (sensitivity,
#' specificity and false positive rate) across a sequence of threshold values.
#' Thresholds are generated from the minimum and maximum predicted values,
#' rounded to the nearest 0.5 by default. For each threshold, predicted values
#' are binarised and compared against observed binary labels using a 2×2
#' confusion matrix. Metrics are returned in a data frame with one row per
#' threshold.
#'
#' @param pred_values Numeric vector of predicted values (e.g. logits or
#' probabilities) used to generate threshold-based classifications.
#' @param obs_values Binary vector of cooresponding observed binary values.
#' Values are coerced to a factor internally (levels = `0` and `1`).
#' @param step Numeric, the increment between threshold values. Defaults to 0.5.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{threshold}{Threshold used to classify predictions.}
#'   \item{sensitivity}{True positive rate (TP / (TP + FN)).}
#'   \item{specificity}{True negative rate (TN / (TN + FP)).}
#'   \item{fpr}{False positive rate (FP / (FP + TN)).}
#' }
#'
#' @examples
#' \dontrun{
#' thresholds_df <- compute_threshold_metrics(
#'   pred_values = predictions$prediction_lodds,
#'   true_labels = predictions$host_status
#' )
#' }
#'
#' @export
compute_threshold_metrics <- function(pred_values,
                                      obs_values,
                                      step = 0.5) {
  # Checks
  if (!is.numeric(pred_values)) {
    stop("'pred_values' must be numeric")
  }
  if (!is.numeric(obs_values) && !is.factor(obs_values)) {
    stop("'observed_values' must be numeric or factor")
  }
  if (length(pred_values) != length(obs_values)) {
    stop("'pred_values' and 'observed_values' must be of the same length")
  }
  if (!is.numeric(step) || length(step) != 1) {
    stop("'step' must be a nueric value of length 1")
  }

  # Ensure binary labels are 0/1
  obv_binary <- factor(obs_values, levels = c(0, 1))

  # Genrate threshold range
  min_val <- floor(min(pred_values) * 2) / 2
  max_val <- ceiling(max(pred_values) * 2) / 2
  thrs <- seq(min_val, max_val, by = step)

  # Initialise vectors
  sensitivity <- numeric(length(thrs))
  specificity <- numeric(length(thrs))
  fpr <- numeric(length(thrs))

  # For every threshold value, calculate sensitivity, specificity and FPR
  for (i in seq_along(thrs)) {
    thr <- thrs[i]

    pred_binary <- factor(ifelse(pred_values > thr, 1, 0), levels = c(0, 1))

    # Confusion matrix
    #            0 (Obs)   1 (Obs)
    # 0 (Pred)   TN        FN
    # 1 (pred)   FP        TP
    conf_matrix <- table(pred_binary, obv_binary)

    # Extract info from matrix
    true_neg <- conf_matrix[1, 1]
    false_pos <- conf_matrix[2, 1]
    false_neg <- conf_matrix[1, 2]
    true_pos <- conf_matrix[2, 2]

    # Compute sensitivity, specificity, and FPR
    sensitivity[i] <- ifelse(true_pos + false_neg == 0, NA,
                             true_pos / (true_pos + false_neg))
    specificity[i] <- ifelse(true_neg + false_pos == 0, NA,
                             true_neg / (true_neg + false_pos))
    fpr[i]  <- ifelse(false_pos + true_neg == 0, NA,
                      false_pos / (false_pos + true_neg))
  }

  # Return data frame
  data.frame(threshold = thrs,
             sensitivity = sensitivity,
             specificity = specificity,
             fpr = fpr)
}


#' Compute the intersection point of sensitivity and specificity curves
#'
#' This function estimates the threshold at which sensitivity and specificity
#' intersect, based on a linear approximation using the two closest points where
#' the curves approach each other. It identifies the two thresholds with the
#' smallest absolute difference between sensitivity and specificity, fits a
#' linear model to each metric across those two points, and solves for their
#' intersection.
#'
#' @param thr_df A \code{data.frame} containing at least the columns
#' \code{sensitivity}, \code{specificity}, and \code{threshold}. The
#' \code{threshold} column should contain the numeric threshold values
#' associated with the sensitivity and specificity estimates.
#'
#' @return A named numeric vector with two elements:
#'   \describe{
#'     \item{x}{The estimated threshold at which sensitivity and specificity
#'              intersect.}
#'     \item{y}{The corresponding value of sensitivity (= specificity) at that
#'              threshold.}
#'   }
#'
#' @examples
#' \dontrun{
#' intercept <- compute_intercept(thresholds_df)
#' intercept
#' }
#'
#' @export
compute_intercept <- function(thr_df) {
  # Checks
  if (!is.data.frame(thr_df)
      || !all(c("sensitivity", "specificity", "threshold")
              %in% colnames(thr_df))) {
    stop("'thresholds_df' must be a data frame containing columns: ",
         "sensitivity, specificity, threshold")
  }

  # Find the two points where lines are closest
  diffs <- abs(thr_df$sensitivity - thr_df$specificity)
  closest <- thr_df[order(diffs)[1:2], ]

  # Fit lines
  lm_sens <- lm(sensitivity ~ threshold, data = closest)
  lm_spec <- lm(specificity ~ threshold, data = closest)

  # Find intercept
  a1 <- unname(coef(lm_sens)["(Intercept)"])
  b1 <- unname(coef(lm_sens)["threshold"])
  a2 <- unname(coef(lm_spec)["(Intercept)"])
  b2 <- unname(coef(lm_spec)["threshold"])

  x_int <- (a2 - a1) / (b1 - b2)
  y_int <- a1 + b1 * x_int
  intercept <- c(threshold = round(x_int, 2), rate = round(y_int, 2))

  intercept
}


#' Compute the true positive rate (TPR)
#'
#' Calculates the true positive rate from a vector of binary predictions and a
#' corresponding vector of binary observations. The true positive rate is
#' defined as the proportion of observed positives correctly predicted.
#'
#' @param predictions Numeric vector of predicted classes, expected to contain
#' only \code{0} and \code{1}.
#' @param observations Numeric vector of observed classes, expected to contain
#' only \code{0} and \code{1}. If not supplied, all observations default to
#' \code{1}. Defaults to \code{rep(1, length(predictions))}.
#'
#' @return A numeric value giving the true positive rate as a percentage,
#' rounded to two decimal places.
#'
#' @details
#' The function assumes that both input vectors contain only binary values,  abd
#' calculates the true positive rate as:
#'
#' \deqn{TPR = \frac{TP}{TP + FN} \times 100}
#'
#' The function assumes that both input vectors contain only binary values.
#'
#' @examples
#' tpr(c(1, 0, 1, 1), c(1, 1, 1, 0))
#'
#' @export
tpr <- function(predictions,
                observations = rep(1, length(predictions))) {
  # Checks
  if (!is.numeric(predictions)) {
    stop("'predictions' must be a numeric vector")
  }
  if (!is.numeric(observations)) {
    stop("'obversations' must be a numeric vector")
  }

  # Confidence matrix
  conf_matrix <- table(
    pred = predictions,
    obs = observations
  )

  # TPR
  true_positives <- conf_matrix["1", "1"]
  false_negatives <- conf_matrix["0", "1"]

  round(true_positives / (true_positives + false_negatives) * 100, 2)
}
