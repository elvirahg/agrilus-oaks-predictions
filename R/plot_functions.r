#' Plot Geographic Distance Metrics
#'
#' Generates one or more scatter plots of distance metrics from a data frame,
#' automatically sorting and transforming them based on their column names.
#' The function supports any combination of raw, normalised, or log-transformed
#' distance columns.
#' Each distance metric is plotted as a scatter plot, and all plots are
#' arranged in a grid layout using \pkg{gridExtra}.
#'
#' @param dist_data A data frame containing the distance metrics to be plotted.
#' @param dist_cols A character vector specifying the column names within
#' `dist_data` to plot. Columns should correspond to distance metrics such as
#'`geo_min_dist_mean`, `geo_min_dist_mean_norm`, or `geo_min_dist_mean_log`.
#'
#' @return
#' A composite grid of \pkg{ggplot2} scatter plots.
#'
#' @examples
#' \dontrun{
#' plot_geo_dists(
#'   interaction_data,
#'   dist_cols = c("geo_min_dist_mean",
#'                 "geo_min_dist_mean_norm",
#'                 "geo_min_dist_mean_log")
#' )
#' }
#'
#' @importFrom gridExtra grid.arrange
#' @export
plot_geo_dists <- function(dist_data, dist_cols) {
  # Checks
  if (!is.data.frame(dist_data)) {
    stop("'dist_data' must be a data frame")
  }
  if (!is.character(dist_cols)) {
    stop("'dist_cols' must be a character vector")
  }
  if (!all(dist_cols %in% colnames(dist_data))) {
    stop("One or more specified columns do not exist in 'dist_data'")
  }

  # Plot colours
  cols <- hcl.colors(length(dist_cols), "Harmonic")

  # Plot list
  dist_plots <- vector("list", length(dist_cols))

  for (i in seq_along(dist_cols)) {
    colname <- dist_cols[i]

    # Plot text
    title_text <- paste0("Distance (", colname, ")")
    if (colname %in% c("geo_min_dist_mean", "geo_min_dist_median")) {
      y_text <- paste0("Distance (km)")
    } else if (colname %in% c("geo_min_dist_mean_norm",
                              "geo_min_dist_median_norm")) {
      y_text <- paste0("Normalised distance")
    } else {
      y_text <- paste0("Log distance")
    }

    # Sort and plot distances
    if (colname %in% c("geo_min_dist_mean", "geo_min_dist_median",
                       "geo_min_dist_mean_log", "geo_min_dist_median_log")) {
      dist_sorted <- sort(dist_data[[colname]], decreasing = FALSE)
      dist_plots[[i]] <- make_dist_plot(dist_sorted, cols[i],
                                        title_text, y_text)
    } else if (colname %in% c("geo_min_dist_mean_norm",
                              "geo_min_dist_median_norm")) {
      dist_sorted <- sort(dist_data[[colname]], decreasing = TRUE)
      dist_plots[[i]] <- make_dist_plot(-dist_sorted, cols[i],
                                        title_text, y_text)
    }
  }

  # Arrange plots
  ncol <- ceiling(length(dist_plots) / 2)
  do.call(gridExtra::grid.arrange, c(dist_plots, ncol = ncol))
}


#' Create a Single Distance Plot
#'
#' Internal helper to generate a scatter plot of a single distance metric.
#'
#' @param values Numeric vector of distances to plot.
#' @param colour Character string specifying point colour.
#' @param title Character string for the plot title.
#' @param ylab_text Character string for the y-axis label.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @import ggplot2
#' @keywords internal
make_dist_plot <- function(values,
                           colour,
                           title,
                           ylab_text) {
  # Checks
  if (!is.numeric(values)) {
    stop("'values' must be a numeric vector")
  }
  if (!is.character(colour)) {
    stop("'colour' must be a character")
  }
  if (!is.character(title)) {
    stop("'title' must be a character")
  }
  if (!is.character(ylab_text)) {
    stop("'ylab_text' must be a character")
  }

  # Plot
  ggplot2::ggplot(data = data.frame(x = seq_along(values), y = values),
                  ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(colour = colour) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Index") +
    ggplot2::ylab(ylab_text)
}


#' Plot comparison of two model prediction columns with a fitted linear
#' relationship
#'
#' This function takes a data frame containing model prediction columns from
#' two different models and a categorical \code{host_status} variable. It fits
#' a simple linear model of \code{full ~ other}, prints a summary, and produces
#' a scatter plot with the fitted line overlaid.
#'
#' @param df A data frame containing the model prediction columns and a
#' \code{host_status} column.
#' @param full A character string giving the name of the prediction column for
#' the “full” model. Must refer to a numeric column in \code{df}.
#' @param other A character string giving the name of the prediction column for
#' the comparison model. Must refer to a numeric column in \code{df}.
#' @param xlab_text A character string used as the x-axis label.
#'
#' @return A \code{ggplot2} object showing the relationship between the two
#' model prediction columns with points coloured by \code{host_status} and a
#' dashed fitted line from the linear model.
#'
#' @details
#' The function is intended for visual comparison of predictions generated by
#' two models. It assumes that both prediction columns are numeric and that
#' \code{host_status} is a factor. The fitted linear model is printed to the
#' console for quick inspection.
#'
#' @examples
#' \dontrun{
#' plot_model_comparison(df = predictions,
#'                       full = "full_model_pred",
#'                       other = "other_model_pred",
#'                       xlab_text = "Other model prediction")
#' }
#'
#' @import ggplot2
#' @export
plot_model_comparison <- function(df,
                                  full,
                                  other,
                                  xlab_text) {
  # Checks
  if (!is.character(full) || length(full) != 1) {
    stop("'full' must be a character of length 1")
  }
  if (!is.character(other) || length(other) != 1) {
    stop("'full' must be a character of length 1")
  }
  if (!is.data.frame(df)
      || !all(c(full, other, "host_status") %in% colnames(df))) {
    stop("'df' must be a data frame containing columns: ",
         paste(full, other, "host_status", collapse = ", "))
  }
  if (!is.numeric(df[[full]]) || !is.numeric(df[[other]])) {
    stop("'df[[full]]' and 'df[[other]]' must be numeric")
  }
  if (!is.factor(df[["host_status"]])) {
    stop("'df[['host_status']]' must be a factor")
  }
  if (!is.character(xlab_text) || (length(xlab_text) != 1)) {
    stop("'xlab_text' must be a character of length 1")
  }

  # Fit the linear model
  lm_formula <- as.formula(paste(full, "~", other))
  lm_summary <- summary(lm(lm_formula, data = df))
  print(lm_summary)

  # Extract intercept and slope
  intercept <- lm_summary$coefficients[1, "Estimate"]
  slope <- lm_summary$coefficients[2, "Estimate"]

  # Plot
  ggplot2::ggplot(df, ggplot2::aes(y = .data[[full]],
                                   x = .data[[other]],
                                   col = host_status)) +
    ggplot2::scale_colour_manual(values = c("lightgrey", "steelblue")) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = intercept,
                         slope = slope,
                         colour = "black",
                         linewidth = 2,
                         linetype = "dashed") +
    ggplot2::ylab("Full model") +
    ggplot2::xlab(xlab_text) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}


#' Plot sensitivity and specificity curves with intersection lines
#'
#' Creates a plot showing sensitivity and specificity as functions of the
#' threshold, and overlays dashed vertical and horizontal lines indicating
#' the intersection point supplied in `thr_ints`.
#'
#' @param thr_df A data frame containing the columns `sensitivity`,
#' `specificity`, and `threshold`. Each row represents one evaluated
#' threshold and its corresponding performance metrics.
#' @param thr_ints A named numeric vector with elements `x` and `y`
#' indicating the threshold (x-coordinate) and rate (y-coordinate) of the
#' intersection point between the sensitivity and specificity lines.
#'
#' @return A ggplot object visualising the sensitivity and specificity
#' curves with the intersection guides.
#'
#' @examples
#' \dontrun{
#'   plot_sens_spec_threshold(thresholds_df, c(x = 0.42, y = 0.87))
#' }
#'
#' @import ggplot2
#' @export
plot_sens_spec_threshold <- function(thr_df, thr_ints) {
  # Checks
  if (!is.data.frame(thr_df)
      || !all(c("sensitivity", "specificity", "threshold")
              %in% colnames(thr_df))) {
    stop("'thresholds_df' must be a data frame containing columns: ",
         "sensitivity, specificity, threshold")
  }
  if (!is.vector(thr_ints)
      || !all(c("threshold", "rate") %in% names(thr_ints))) {
    stop("'thr_intercepts' must be a named vector containing",
         "'threshold' and 'rate'")
  }
  # Plot
  ggplot2::ggplot(thr_df, ggplot2::aes(threshold)) +
    ggplot2::geom_line(aes(y = sensitivity),
                       colour = "steelblue", linewidth = 2) +
    ggplot2::geom_line(aes(y = specificity),
                       colour = "coral", linewidth = 2) +
    ggplot2::geom_vline(xintercept = thr_ints["threshold"],
                        linetype = "dashed", colour = "darkgrey") +
    ggplot2::geom_hline(yintercept = thr_ints["rate"],
                        linetype = "dashed", colour = "darkgrey") +
    ggplot2::labs(x = "Threshold", y = "Rate") +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot predicted log-odds by host status with violin and box plots
#'
#' Produces a violin plot (with an embedded box plot) of predicted log-odds
#' for each host status category, with an optional horizontal line marking a
#' chosen threshold, if probided.
#'
#' @param predictions A data frame containing at least two columns:
#'   \itemize{
#'     \item{prediction_lodds}{Numeric vector of predicted log-odds.}
#'     \item{host_status}{Factor or character vector indicating class labels
#'       (e.g., 0/1 or negative/positive).}
#'   }
#'
#' @param threshold A single numeric value giving the log-odds threshold to draw
#' as a horizontal line, or \code{NULL} (default) to omit the line.
#'
#' @return A \code{ggplot} object representing the violin/box plot with a
#' threshold line.
#'
#' @examples
#' \dontrun{
#' plot_predictions_violin(predictions_df, threshold = 0.5)
#' }
#'
#' @import ggplot2
#' @export
plot_predictions_violin <- function(predictions,
                                    threshold = NULL) {
  # Checks
  if (!is.data.frame(predictions)
      || !all(c("prediction_lodds", "host_status")
              %in% colnames(predictions))) {
    stop("'predictions' must be a data frame containing 'prediction_lodds'",
         "and 'host_status'")
  }
  if (!is.null(threshold)
      && (!is.numeric(threshold) || length(threshold) != 1)) {
    stop("'threshold' must be a single numeric value or NULL")
  }
  # Make sure host status is a factor
  predictions$host_status <- as.factor(predictions$host_status)

  # Plot
  p <- ggplot2::ggplot(data = predictions,
                       ggplot2::aes(y = prediction_lodds,
                                    x = host_status,
                                    fill = host_status)) +
    ggplot2::geom_violin(width = 1) +
    ggplot2::geom_boxplot(width = 0.1, fill = "grey") +
    ggplot2::scale_fill_manual(values = c("aquamarine2", "coral")) +
    ggplot2::ggtitle("Predicted interactions") +
    ggplot2::xlab("Host status") +
    ggplot2::ylab("Log-odds predicted host status") +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 16),
                   axis.title = ggplot2::element_text(size = 16))

  # Add threshold if provided
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold,
                                 linewidth = 1.5,
                                 linetype = "dotted",
                                 col = "maroon")
  }

  p
}


#' Plot a fourfold confusion matrix and optionally print performance metrics
#'
#' Generates a fourfold plot of a binary confusion matrix using predicted and
#' observed class labels. The function also computes and prints sensitivity
#' (true positive rate) and specificity (true negative rate) if requested.
#' Input vectors may be factors or numeric, containing only 0s and 1s.
#'
#' @param observed_status A factor or numeric vector containing observations,
#' coded as 0 and 1.
#' @param prediction_binary A factor or numeric vector containing predictions
#' , coded as 0 and 1.
#' @param print_metrics Logical; if \code{TRUE} (default), sensitivity and
#' specificity are calculated and printed.
#'
#' @return Generates a fourfold plot and, optionally, prints classification
#' metrics.
#'
#' @examples
#' \dontrun{
#' obs <- c(0,1,1,0,1)
#' pred <- c(0,1,0,0,1)
#' plot_fourfold(obs, pred)
#' }
#'
#' @export
plot_fourfold <- function(observations,
                          predictions,
                          print_metrics = TRUE) {
  if (!is.factor(observations)
      && !is.numeric(observations)
      && !is.integer(observations)) {
    stop("'observations' must be a factor or numeric vector containing",
         "0 and 1")
  }
  if (!all(observations %in% c("0", "1", 0, 1))) {
    stop("'observations' must contain only values 0 and 1")
  }
  if (!is.factor(predictions)
      && !is.numeric(predictions)
      && !is.integer(predictions)) {
    stop("'prediction_binary' must be a factor or numeric vector containing",
         "0 and 1")
  }
  if (!all(predictions %in% c("0", "1", 0, 1))) {
    stop("'prediction_binary' must contain only values 0 and 1")
  }

  # Confusion matrix
  conf_matrix <- table(prediction = predictions,
                       observed   = observations)

  if (print_metrics) {
    # Extract metrics
    true_neg <- conf_matrix["0", "0"]
    false_pos <- conf_matrix["1", "0"]
    false_neg <- conf_matrix["0", "1"]
    true_pos <- conf_matrix["1", "1"]

    tp_rate <- round(true_pos / (true_pos + false_neg) * 100, 2)
    tn_rate <- round(true_neg / (true_neg + false_pos) * 100, 2)

    # Print results
    cat(paste("True positive rate (%):", tp_rate, "\n"))
    cat(paste("True negative rate (%):", tn_rate, "\n"))
  }

  # Fourfold plot
  fourfoldplot(conf_matrix,
               color = c("#CC6666", "#99CC99"),
               conf.level = 0,
               main = "Confusion Matrix",
               margin = c(2))
}


#' Plot predicted values grouped by species
#'
#' Creates a faceted scatter plot of predicted values (e.g. log-odds)
#' against host status, grouped by a species column. An optional threshold
#' line can be added, and the plot may be restricted to known host species
#' (i.e., `host_status == 1`).
#'
#' @param data A data frame containing at least `host_status`,
#' the species column specified in `species_col`, and the variable named in
#' `pred_var`.
#' @param species_col A single character string giving the name of the column
#' defining species groups for faceting.
#' @param threshold Optional single numeric value. If supplied, a horizontal
#' line is drawn at this value. Default is `NULL`.
#' @param pred_var A single character string giving the name of the numeric
#' variable to plot on the y-axis. Default is `"prediction_lodds"`.
#' @param hosts_only Logical. If `TRUE`, the function restricts the data to
#' species known host species (`host_status == 1`). Default is `FALSE`.
#' @param ncol Single integer giving the number of facet columns. Default is 16.
#' @param point_size Single numeric value indicating the size of plotted points.
#' Default is 1.
#' @param title Optional single character string for the plot title. Default is
#' `NULL`.
#' @param italic_text Logical indicating whether species labels in facet strips
#' are printed in italic. Default is `TRUE`.
#'
#' @return A `ggplot` object.
#'
#' @details
#' The plot displays host status on the x-axis and the chosen predictive
#' variable on the y-axis, faceted by species.
#'
#' @examples
#' \dontrun{
#' plot_predictions_by_species(
#'   data = predictions,
#'   species_col = "quercus_sp",
#'   threshold = 0.5,
#'   ncol = 12,
#'   title = "Predicted interactions by species"
#' )
#' }
#'
#' @import ggplot2
#' @export
plot_predictions_by_species <- function(data,
                                        species_col,
                                        threshold = NULL,
                                        pred_var = "prediction_lodds",
                                        hosts_only = FALSE,
                                        ncol = 16,
                                        point_size = 1,
                                        title = NULL,
                                        italic_text = TRUE) {
  # Checks
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!all(c("host_status", species_col) %in% names(data))) {
    stop(paste("'data' must contain host_status and", species_col))
  }
  if (!is.character(species_col) || length(species_col) != 1) {
    stop("'species_col' must be a single character string")
  }
  if (!is.null(threshold)
      && (!is.numeric(threshold) || length(threshold) != 1)) {
    stop("'threshold' must be a single numeric value or NULL")
  }
  if (!is.character(pred_var) || length(pred_var) != 1) {
    stop("'y_var' must be a single character string")
  }
  if (!is.logical(hosts_only) || length(hosts_only) != 1) {
    stop("'hosts_only' must be a single logical value")
  }
  if (!is.numeric(ncol) || length(ncol) != 1 || ncol != as.integer(ncol)) {
    stop("'ncol' must be a single integer")
  }
  if (!is.numeric(point_size) || length(point_size) != 1) {
    stop("'point_size' must be a single numeric value")
  }
  if (!is.null(title)
      && (!is.character(title) || length(title) != 1)) {
    stop("'title' must be a single character string or NULL")
  }
  if (!is.logical(italic_text) || length(italic_text) != 1) {
    stop("'italic_text' must be a single logical value")
  }

  # If host_only, only retain host species
  if (hosts_only) {
    host_species <- unique(data[data$host_status == 1, ][[species_col]])
    data <- data[data[[species_col]] %in% host_species, ]
  }

  # Plot
  species_formula <- as.formula(paste("~", species_col))

  p <- ggplot2::ggplot(data,
                       ggplot2::aes(x = host_status,
                                    y = .data[[pred_var]])) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(species_formula, ncol = ncol)

  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold,
                                 linewidth = 1.2,
                                 linetype = "dotted",
                                 colour = "maroon")
  }

  if (italic_text) {
    p <- p + ggplot2::theme(strip.text = ggplot2::element_text(face = "italic"))
  }

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  p
}


#' Plot predicted host interactions and highlight species above a threshold
#'
#' This function visualises predicted interactions for a single hosted species
#' by plotting prediction scores against host status.
#' Host-taxon species whose predicted values exceed a given threshold (but are
#' observed non-hosts) are highlighted with text labels.
#'
#' @param data A data frame containing prediction results.
#' @param hosted_col Character string specifying the column name for the hosted
#' species.
#' @param host_col Character string specifying the column name for the host-
#' taxon species.
#' @param pred_col Character string giving the column name of the prediction
#' scores (default: `"prediction_lodds"`).
#' @param threshold A single numeric value; predictions equal to or above this
#' are considered to exceed the threshold.
#' @param size A single integer giving the point size for plotted points.
#' Default is 5.
#' @param label_size A single integer giving the text size for highlighted
#' labels. Default is 4.
#'
#' @return A ggplot object showing prediction scores, with species above the
#' threshold highlighted.
#'
#' @examples
#' \dontrun{
#' plot_highlighted_predictions(
#'   data = predictions,
#'   hosted_col = "agrilus_sp",
#'   host_col = "quercus_sp",
#'   threshold = 0.5
#' )
#' }
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @export
plot_highlighted_predictions <- function(data,
                                         hosted_col,
                                         host_col,
                                         pred_col = "prediction_lodds",
                                         threshold,
                                         size = 5,
                                         label_size = 4) {
  # Checks
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!all(c("host_status", hosted_col, host_col, pred_col) %in% names(data))) {
    stop(paste("'data' must contain host_status", hosted_col,
               host_col, pred_col, sep = ", "))
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("'threshold' must be a single numeric value.")
  }
  if (!is.numeric(size) || length(size) != 1 || size != as.integer(size)) {
    stop("'size' must be a single integer")
  }
  if (!is.numeric(label_size) || length(label_size) != 1
      || label_size != as.integer(label_size)) {
    stop("'label_size' must be a single integer")
  }

  # Find predicted hosts for focal hosted species, i.e., non-host, host taxon
  # species above threshold (where host_status == 0)
  pred_hosts <- data[data[["prediction_lodds"]] >= threshold, ]
  pred_hosts <- pred_hosts[pred_hosts$host_status == 0, ][[host_col]]

  # Plot
  ggplot2::ggplot(data,
                  ggplot2::aes(x = host_status,
                               y = .data[[pred_col]])) +
    ggplot2::geom_point(size = size) +
    ggplot2::facet_wrap(as.formula(paste("~", hosted_col)), ncol = 8) +
    ggplot2::geom_hline(yintercept = threshold,
                        linetype = "dashed",
                        colour = "red") +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = ifelse(.data[[host_col]] %in% pred_hosts,
                                  .data[[host_col]], "")),
      hjust = 0,
      vjust = 0,
      max.overlaps = 100,
      size = label_size
    )
}


#' Plot changes in predicted values relative to a threshold
#'
#' Produces a scatter plot comparing original and updated prediction values,
#' with colour highlighting based either on whether predictions cross a
#' threshold or on supplied observation values. A prediction is considered to
#' have changed if it is greater than or equal to \code{threshold} in
#' \code{pred_original} and below \code{threshold} in \code{pred_new}.
#' Optionally prints the rows corresponding to such changes.
#'
#' @param pred_original Numeric vector of original prediction values.
#' @param pred_new Numeric vector of updated prediction values.
#' @param threshold Single numeric threshold used to identify prediction
#' changes.
#' @param obs_vals Optional numeric vector of observation values used for
#' colouring when \code{col_by = "obs_vals"}. Must be the same length as
#' \code{pred_original}.
#' @param print_change Logical; if \code{TRUE}, prints rows where predictions
#' cross the threshold from above to below. Defaults to \code{TRUE}.
#' @param col_by Character string indicating the variable used for colour:
#' \code{"pred_change"} (default) colours points by threshold crossing,
#' whereas \code{"obs_vals"} colours points using \code{obs_vals}.
#' @param x_lab Character string for the x-axis label. Defaults to
#' \code{"Original predictions"}.
#' @param y_lab Character string for the y-axis label. Defaults to
#' \code{"New predictions"}.
#' @param title Character string for the plot title. Defaults to
#' \code{"Prediction change"}.
#'
#' @return A \pkg{ggplot2} object showing original vs. new predictions.
#'
#' @details
#' When \code{col_by = "pred_change"}, points are coloured according to whether
#' they cross the threshold. When \code{col_by = "obs_vals"}, the supplied
#' observation values determine the colouring instead. The diagonal reference
#' line is added to help visualise deviations between original and new
#' predictions.
#'
#' @examples
#' \dontrun{
#' plot_prediction_change(
#'   pred_original = c(-3, -1, 0, 2),
#'   pred_new = c(-4, -2, -0.5, 1.5),
#'   threshold = 0,
#'   col_by = "pred_change"
#' )
#' }
#'
#' @import ggplot2
#' @export
plot_prediction_change <- function(pred_original,
                                   pred_new,
                                   threshold,
                                   obs_vals = NULL,
                                   print_change = TRUE,
                                   col_by = "pred_change",
                                   x_lab = "Original predictions",
                                   y_lab = "New predictions",
                                   title = "Prediction change") {
  # Checks
  if (!is.numeric(pred_original)) {
    stop("'pred_original' must be a numeric vector")
  }
  if (!is.numeric(pred_new)) {
    stop("'pred_new' must be a numeric vector")
  }
  if (length(pred_original) != length(pred_new)) {
    stop("'pred_original' and 'pred_new' must be of equal length")
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("'threshold' must be a single numeric value")
  }
  if (!is.logical(print_change) || length(print_change) != 1) {
    stop("'print_change' must be a single logical value")
  }
  if (col_by != "pred_change" && col_by != "obs_vals") {
    stop("'col_by' must be one of 'pred_change' or 'obs_vals'")
  }
  if (col_by == "obs_vals" && is.null(obs_vals)) {
    stop("if 'col_by = 'obs_vals'', obs_val must be a numeric vector")
  }
  if (!(is.null(obs_vals))) {
    if (col_by != "obs_vals") {
      warning("'obs_vals' object will not be used, as col_by != 'obs_vals'")
    }
    if (!(is.numeric(obs_vals)) || length(obs_vals) != length(pred_original)) {
      stop("'obs_vals' must be a numeric vector of equal length to",
           "'pred_original'")
    }
  }
  if (!is.character(x_lab) || length(x_lab) != 1) {
    stop("'x_lab' must be a single character string")
  }
  if (!is.character(y_lab) || length(y_lab) != 1) {
    stop("'y_lab' must be a single character string")
  }
  if (!is.character(title) || length(title) != 1) {
    stop("'title' must be a single character string")
  }

  # Construct a data frame for plotting and printing
  # pred_change: TRUE if original prediction > threshold but new < threshold
  preds_df <- data.frame(
    pred_original = pred_original,
    pred_new = pred_new,
    pred_change = pred_original >= threshold & pred_new < threshold
  )
  if (col_by == "obs_vals") {
    preds_df$obs_vals <- as.logical(obs_vals)
  }

  # Print change
  if (print_change) {
    if (any(preds_df$pred_change)) {
      print(preds_df[preds_df$pred_change, ])
    } else {
      cat("No predictions moved across threshold\n")
    }
  }

  # Plot
  if (col_by == "pred_change") {
    p <- ggplot2::ggplot(preds_df,
                         ggplot2::aes(x = pred_original,
                                      y = pred_new,
                                      col = pred_change))
  } else if (col_by == "obs_vals") {
    p <- ggplot2::ggplot(preds_df,
                         ggplot2::aes(x = pred_original,
                                      y = pred_new,
                                      col = obs_vals))
  }
  p <- p +
    ggplot2::geom_abline(linetype = "dashed",
                         col = "darkgrey",
                         linewidth = 1.5) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("black", "coral")) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.position = "none")

  p
}


#' Plot phylogeny with predicted and known hosts
#'
#' This function takes a phylogenetic tree (`phylo`) and a data frame containing
#' information on predicted and known hosts for each tip, and produces a
#' circular (or rectangular) tree with:
#' - Predicted hosts highlighted with bold italic labels and custom colour.
#' - Known hosts marked with a tip point and custom colour.
#'
#' @param phylo An object of class `phylo` representing the tree to plot.
#' @param pred_hosts A data.frame containing tip information with three columns:
#'   - `label`: character vector of tip names (species names matching tip
#' labels).
#'   - `is_known`: logical vector indicating if the tip is a known host.
#'   - `is_pred`: logical vector indicating if the tip is a predicted host.
#' @param label Character string indicating the column name in `pred_hosts`
#'   containing the tip labels.
#' @param layout Character string specifying the tree layout (default:
#' `circular`).
#' @param tiplab_offset Numeric, offset distance of tip labels from the tree
#' (defaut: `5`).
#' @param tiplab_size Numeric, font size of tip labels (default: 1.7).
#' @param tippoint_shape Numeric, shape code for known host tip points (default:
#' `18`).
#' @param tippoint_size Numeric, size of the tip point symbols (default: 2).
#' @param pred_color Character, colour for predicted host labels (default:
#' `darkgoldenrod3`).
#' @param nonpred_color Character, colour for non-predicted tips (default:
#' `gray50`).
#' @param known_color Character, colour for known host tip points (default:
#' `#cf3530`).
#' @param show_legend Logical, whether to show a legend (default: `FALSE`).
#'
#' @return A `ggtree` object that can be printed or further customised.
#'
#' @examples
#' \dontrun{
#' pred_hosts <- data.frame(
#'   quercus_sp = oak_phylo$tip.label,
#'   is_known = oak_phylo$tip.label %in% known_hosts,
#'   is_pred = oak_phylo$tip.label %in% pred_hosts
#' )
#' plot_oak_phylo_hosts(oak_phylo, pred_hosts, label = "quercus_sp")
#' }
#'
#' @import ggplot2 ggtree treeio
#' @export
plot_phylo_hosts <- function(phylo,
                             pred_hosts,
                             label,
                             layout = "circular",
                             tiplab_offset = 5,
                             tiplab_size = 1.7,
                             tippoint_shape = 18,
                             tippoint_size = 2,
                             pred_color = "darkgoldenrod3",
                             nonpred_color = "gray50",
                             known_color = "#cf3530",
                             show_legend = FALSE) {
  # Checks
  if (!is.data.frame(pred_hosts)
      || !(all(c(label, "is_known", "is_pred") %in% colnames(pred_hosts)))) {
    stop(paste0("'pred_hosts' must be a data.frame with colnames: ",
                label, ", is_known, is_pred"))
  }
  if (!is.character(pred_hosts[[label]])) {
    stop("'pred_hosts[[label]]' must be a character vector")
  }
  if (!is.logical(pred_hosts$is_known)) {
    stop("'pred_hosts$is_known' must be a logical vector")
  }
  if (!is.logical(pred_hosts$is_pred)) {
    stop("'pred_hosts$is_pred' must be a logical vector")
  }
  if (class(phylo) != "phylo") {
    stop("'phylo' must be an object of class 'phylo'")
  }
  if (!is.character(label) || length(label) != 1) {
    stop("'label' must be a character string of length 1")
  }
  if (!is.character(layout) || length(layout) != 1) {
    stop("'layout' must be a character string of length 1")
  }
  if (!is.numeric(tiplab_offset) || length(tiplab_offset) != 1) {
    stop("'tiplab_offset' must be a single numeric value")
  }
  if (!is.numeric(tiplab_size) || length(tiplab_size) != 1) {
    stop("'tiplab_size' must be a single numeric value")
  }
  if (!is.numeric(tippoint_shape) || length(tippoint_shape) != 1) {
    stop("'tippoint_shape' must be a single numeric value")
  }
  if (!is.numeric(tippoint_size) || length(tippoint_size) != 1) {
    stop("'tippoint_size' must be a single numeric value")
  }
  if (!is.character(pred_color) || length(pred_color) != 1) {
    stop("'pred_color' must be a character string of length 1")
  }
  if (!is.character(nonpred_color) || length(nonpred_color) != 1) {
    stop("'nonpred_color' must be a character string of length 1")
  }
  if (!is.character(known_color) || length(known_color) != 1) {
    stop("'known_color' must be a character string of length 1")
  }
  if (!is.logical(show_legend) || length(show_legend) != 1) {
    stop("'show_legend' must be a single logical value")
  }

  # Add host info to phylo
  phylo_preds <- phylo |>
    ggtree::fortify() |>
    dplyr::left_join(pred_hosts, by = c("label" = label)) |>
    treeio::as.treedata()

  # Plot
  ggtree::ggtree(phylo_preds, layout = layout) +
    # Tip labels for predicted hosts
    ggtree::geom_tiplab(ggplot2::aes(subset = TRUE,
                                     fontface = ifelse(is_pred,
                                                       "bold.italic",
                                                       "italic"),
                                     color = is_pred),
                        offset = tiplab_offset,
                        size = tiplab_size,
                        show.legend = show_legend) +
    ggplot2::scale_color_manual(values = c(nonpred_color, pred_color)) +
    # Tip points for known hosts
    ggtree::geom_tippoint(ggplot2::aes(subset = is_known),
                          shape = tippoint_shape,
                          size = tippoint_size,
                          color = known_color,
                          show.legend = show_legend)
}


#' Plot plant species distribution per region (level)
#'
#' This function prepares counts of plant species per region and plots them
#' on a world map with discrete colour bins.
#'
#' @param df A data frame (or sf object) containing plant distribution data.
#' @param group_level Character string specifying the column to group by
#' (default `"LEVEL3_NAM"`).
#' @param counts_type Character string specifying which counts to plot.
#'   Options are `"sp"`, `"native_sp"`, `"hosts"`, `"hosts_native"`, `"pred"`,
#' `"native_pred"`. Defaults to `"species"`.
#' @param world_sf Optional `sf` object with world polygons. If `NULL`
#' (default), uses rnaturalearth.
#' @param breaks Numeric vector specifying breakpoints for binning counts.
#' @param labels Character vector of labels for the bins.
#' @param palette_val Colour palette values for `scale_fill_brewer`; must be a
#' character vector. Defaults to RColorBrewer::brewer.pal(9, "PuBuGn").
#' @param na_col Colour for NA values. Defaults to `"gray70"`.
#'
#' @return A ggplot object showing a world map coloured by the selected counts.
#'
#' @seealso \code{\link{prepare_counts_region}}
#'
#' @examples
#' \dontrun{
#' plot_distribution_region(oak_distributions,
#'                          counts_type = "native_sp")
#' }
#'
#' @import rnaturalearth ggplot2 RColorBrewer
#' @export
plot_distribution_region <- function(df,
                                     group_level = "LEVEL3_NAM",
                                     counts_type = "species",
                                     world_sf = NULL,
                                     breaks = c(-Inf, 0, 10, 20,
                                                30, 40, Inf),
                                     labels = c("0", "1–10", "11–20",
                                                "21–30", "31–40", "40+"),
                                     palette_vals = RColorBrewer::brewer.pal(9, "PuBuGn"),
                                     na_col = "gray70") {
  # Checks
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.null(world_sf) && !inherits(world_sf, "sf")) {
    stop("'world_sf' must be a world sf object")
  }
  if (!is.character(palette_vals)) {
    stop("'palette' must be a character vector")
  }


  # Prepare world map if not provided
  if (is.null(world_sf)) {
    world_sf <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf")
  }

  # Prepare counts (with binned values)
  counts_binned <- prepare_counts_region(df = df,
                                         group_level = group_level,
                                         counts_type = counts_type,
                                         breaks = breaks,
                                         labels = labels)

  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = world_sf,
                     fill = na_col,
                     colour = "gray40",
                     size = 0.2) +
    ggplot2::geom_sf(data = counts_binned,
                     ggplot2::aes(fill = value_band),
                     colour = "gray40",
                     size = 0.1) +
    ggplot2::scale_fill_manual(values = palette_vals,
                               na.value = na_col) +
    ggplot2::theme_minimal() +
    ggplot2::labs(fill = counts_type,
                  title = paste("Choropleth map of", counts_type,
                                "distribution, grouped by", group_level))

}


#' Prepare counts per region (level) and bin them
#'
#' This is an internal helper function that computes the number of
#' plant species, native species, known hosts, or predicted hosts
#' per region, and assigns each region to a discrete bin for plotting.
#'
#' @param df A data frame (or sf object) containing plant distribution
#' data. Must include the columns `plant_sp`, `occurrence_type`,
#' `known_host` (for host couts), `pred_host` (for predicte host counts),
#' and the grouping column.
#' @param group_level Character string giving the name of the column to group
#' by. Defaults to `"LEVEL3_NAM"`.
#' @param counts_type Character string specifying what to count. Options are:
#' `"species"`, `"species_native"`, `"hosts"`, `"hosts_native"`, `"pred_hosts"`,
#' `"pred_hosts_native"`.
#' @param breaks Numeric vector of break points for binning counts. Defaults to
#' `c(-Inf, 0, 10, 20, 30, 40, Inf)`.
#' @param labels Character vector of labels for the bins; must have length one
#' less than `breaks`. Defaults to
#' `c("0", "1–10", "11–20", "21–30", "31–40", "40+")`.
#'
#' @return A data frame with columns:
#'   - the grouping column (`group_level`)
#'   - `value`: the raw count for the selected `counts_type`
#'   - `value_band`: the binned factor for plotting
#'
#' @import dplyr
#' @keywords internal
prepare_counts_region <- function(df,
                                  group_level = "LEVEL3_NAM",
                                  counts_type,
                                  breaks = c(-Inf, 0, 10, 20, 30, 40, Inf),
                                  labels = c("0", "1–10", "11–20", "21–30",
                                             "31–40", "40+")) {
  # Checks
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.character(group_level) || length(group_level) != 1) {
    stop("'group_level' must be a single character string")
  }
  if (!is.character(counts_type) || length(counts_type) != 1) {
    stop("'counts_type' must be a single character string")
  }
  if (!(counts_type %in% c("species", "species_native",
                           "hosts", "hosts_native",
                           "pred_hosts", "pred_hosts_native"))) {
    stop("'counts_type' must be one of 'species', 'species_native',
         'hosts', 'hosts_native', 'pred_hosts', or 'pred_hosts_native'")
  }
  if (!is.numeric(breaks)) {
    stop("'breaks' must be a numeric vector")
  }
  if (!is.character(labels)) {
    stop("'labels' must be a character vector")
  }

  # Generate counts data frame
  counts <- df |>
    dplyr::group_by(!!dplyr::sym(group_level)) |>
    dplyr::summarise(
      value = dplyr::case_when(
        counts_type == "species" ~
          n_distinct(plant_sp),
        counts_type == "species_native" ~
          n_distinct(plant_sp[occurrence_type == "native"]),
        counts_type == "hosts" ~
          n_distinct(plant_sp[known_host]),
        counts_type == "hosts_native" ~
          n_distinct(plant_sp[known_host & occurrence_type == "native"]),
        counts_type == "pred_hosts" ~
          n_distinct(plant_sp[pred_host]),
        counts_type == "pred_hosts_native" ~
          n_distinct(plant_sp[pred_host & occurrence_type == "native"])
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(!!dplyr::sym(group_level)) |>
    dplyr::mutate(value_band = cut(value,
                                   breaks = breaks,
                                   labels = labels))

  cat("Counts value summary:\n")
  summary_result <- summary(counts$value)
  cat(paste(names(summary_result), summary_result, sep = ": "), sep = "\n")

  counts
}


#' Plot plant species distribution per country
#'
#' This function computes counts of plant species per country and plots
#' them on a world map with discrete colour bins.
#'
#' @param df An sf object or data frame containing plant distribution data.
#' Must include the columns `plant_sp`, `occurrence_type`, `known_host` (for
#' host counts), `pred_host` (for predicted host counts).
#' @param counts_type Character string specifying which counts to plot.
#' Options are `"species"`, `"species_native"`, `"hosts"`, `"hosts_native"`,
#' `"pred_hosts"`, `"pred_hosts_native"`. Defaults to `"species"`.
#' @param world_sf Optional `sf` object with world country polygons. If `NULL`
#' (default), the world map is obtained from `rnaturalearth::ne_countries()`.
#' @param country_lookup Optional data frame used to map distribution regions
#' to country names. Must contain columns `LEVEL3_COD` and `COUNTRY`.
#' @param breaks Numeric vector specifying breakpoints for binning counts.
#' Defaults to `c(1, seq(10, 90, 10))`.
#' @param labels Character vector of labels for the bins. Defaults to
#' `c(paste0(seq(1, 90, 10), "–", seq(10, 90, 10)))`.
#' @param palette_vals Colour palette values for `scale_fill_manual`; must be
#' a character vector. Defaults to
#' `colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(9)`.
#' @param na_col Colour for NA values. Defaults to `"gray70"`.
#'
#' @return A ggplot object showing a world map coloured by the selected counts.
#'
#' @seealso \code{\link{prepare_counts_country}}
#'
#' @examples
#' \dontrun{
#' plot_distribution_country(oak_distributions,
#'                           counts_type = "species_native")
#' }
#'
#' @import ggplot2 RColorBrewer
#' @export
plot_distribution_country <- function(df,
                                      counts_type = "species",
                                      world_sf = NULL,
                                      country_lookup = NULL,
                                      breaks = c(1, seq(10, 90, 10)),
                                      labels = c(paste0(seq(1, 90, 10),
                                                        "–",
                                                        seq(10, 90, 10))),
                                      palette_vals = colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(9),
                                      na_col = "gray70") {
  # Checks
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.character(palette_vals)) {
    stop("'palette' must be a character vector")
  }

  # Get world counts
  world_counts <- prepare_counts_country(df = df,
                                         counts_type = counts_type,
                                         world_sf = world_sf,
                                         country_lookup = country_lookup,
                                         breaks = breaks,
                                         labels = labels)

  ggplot2::ggplot(world_counts) +
    ggplot2::geom_sf(ggplot2::aes(fill = value_band)) +
    ggplot2::scale_fill_manual(values = palette_vals,
                               na.value = na_col) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      fill = counts_type,
      title = paste("Country-level choropleth of", counts_type)
    )
}


#' Prepare counts per country
#'
#' This is an internal helper function that computes the number of plant
#' species, native species, known hosts, or predicted hosts per country
#' and assigns each country to a discrete bin for plotting.
#'
#' @param df A data frame (or sf object) containing plant distribution
#' data. Must include the columns `plant_sp`, `occurrence_type`,
#' `known_host` (for host counts), `pred_host` (for predicted host counts).
#' @param counts_type Character string specifying what to count. Options are:
#' `"species"`, `"species_native"`, `"hosts"`, `"hosts_native"`,
#' `"pred_hosts"`, `"pred_hosts_native"`. Defaults to `"species"`.
#' @param world_sf Optional `sf` object with world country polygons. If `NULL`,
#' the world map is obtained from `rnaturalearth::ne_countries()`.
#' @param country_lookup Optional data frame used to map distribution regions
#' to country names. Must contain columns `LEVEL3_COD` and `COUNTRY`. No
#' default, but rWCVP::wgsrpd_mapping can be used to generate lookup data frame.
#' @param world_col_map Character string giving the column name in `world_sf`
#' used for joining country names (e.g. `"subunit"`). Defaults to `"subunit"`.
#' @param breaks Numeric vector of break points for binning counts. Defaults to
#' `c(0, 1, seq(10, 90, 10))`.
#' @param labels Character vector of labels for the bins; must have length one
#' less than `breaks`. Defaults to
#' `c(paste0(seq(1, 90, 10), "–", seq(10, 90, 10)))`.
#'
#' @return An sf object of countries with two additional columns:
#'   - `species_count`: raw count for the selected `counts_type`
#'   - `value_band`: the binned factor for plotting
#'
#' @import sf dplyr rnaturalearth
#' @keywords internal
prepare_counts_country <- function(df,
                                   counts_type = "species",
                                   world_sf = NULL,
                                   country_lookup = NULL,
                                   world_col_map = "subunit",
                                   breaks = c(1, seq(10, 90, 10)),
                                   labels = c(paste0(seq(1, 90, 10),
                                                     "–",
                                                     seq(10, 90, 10)))) {
  # Checks
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.character(counts_type) || length(counts_type) != 1) {
    stop("'counts_type' must be a single character string")
  }
  if (!(counts_type %in% c("species", "species_native",
                           "hosts", "hosts_native",
                           "pred_hosts", "pred_hosts_native"))) {
    stop("'counts_type' must be one of 'species', 'species_native',
         'hosts', 'hosts_native', 'pred_hosts', or 'pred_hosts_native'")
  }
  if (!is.null(world_sf) && !inherits(world_sf, "sf")) {
    stop("'world_sf' must be a world sf object")
  }
  if (!is.null(country_lookup) && !inherits(country_lookup, "data.frame")) {
    stop("'country_lookup' must be a data frame")
  }
  if (!(any(c("LEVEL3_COD", "COUNTRY") %in% colnames(country_lookup)))) {
    stop("'country_lookup' must contain columns 'LEVEL3_COD' and 'COUNTRY'")
  }
  if (!is.numeric(breaks)) {
    stop("'breaks' must be a numeric vector")
  }
  if (!is.character(labels)) {
    stop("'labels' must be a character vector")
  }

  # Prepare world map if not provided
  if (is.null(world_sf)) {
    world_sf <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf")
  }

  # Add country info to df
  df_country <- df |>
    dplyr::left_join(country_lookup,
                     by = "LEVEL3_COD",
                     relationship = "many-to-many")

  # Compute counts according to type
  country_counts <- df_country |>
    dplyr::filter(!is.na(COUNTRY)) |>
    dplyr::group_by(COUNTRY) |>
    dplyr::summarise(
      species_count = dplyr::case_when(
        counts_type == "species" ~
          dplyr::n_distinct(plant_sp),

        counts_type == "species_native" ~
          dplyr::n_distinct(plant_sp[occurrence_type == "native"]),

        counts_type == "hosts" ~
          dplyr::n_distinct(plant_sp[known_host]),

        counts_type == "hosts_native" ~
          dplyr::n_distinct(plant_sp[known_host & occurrence_type == "native"]),

        counts_type == "pred_hosts" ~
          dplyr::n_distinct(plant_sp[pred_host]),

        counts_type == "pred_hosts_native" ~
          dplyr::n_distinct(
            plant_sp[pred_host & occurrence_type == "native"]
          )
      ),
      .groups = "drop"
    )

  # Add country counts to world map
  world_sf <- world_sf |>
    dplyr::left_join(sf::st_drop_geometry(country_counts),
                     by = stats::setNames("COUNTRY", world_col_map))

  # Add value bands to world map
  world_sf$value_band <- cut(world_sf$species_count,
                             breaks = breaks,
                             labels = labels,
                             right = FALSE)

  cat("Counts value summary:\n")
  summary_result <- summary(world_sf$species_count[!is.na(world_sf$species_count)])
  cat(paste(names(summary_result), summary_result, sep = ": "), sep = "\n")

  world_sf
}


#' Plot a phylogeny highlighting host species for a focal species
#'
#' This function plots a phylogenetic tree and highlights host species
#' associated with a given focal species using tip label colours.
#'
#' @param phylo An object of class \code{phylo} containing the phylogenetic
#' tree.
#' @param host_spp A character string specifying the host species. This can be
#' provided as an alternative to inferring the host species via `get_host_spp`
#' when an `interaction_df` is supplied.  Default: `NULL` (not required if
#' `interaction_df` is provided).
#' @param interaction_df A data frame containing species-host interaction data.
#' Default: `NULL` (not required if `host_spp` is provided).
#' @param species Character string giving the focal species of interest.
#' @param species_col A character string indicating the column name in
#' `interaction_df` that contains species names.  Default: `NULL` (not required
#' if `interaction_df` is provided).
#' @param host_col A character string indicating the column name in both `model`
#' and `interaction_df` that contains host species names.
#' @param host_status_col A character string specifying the column in
#' `interaction_df` that indicates the host status.  Default: `NULL` (not
#' required if `interaction_df` is provided).
#' @param main Character string indicating plot title. Defaults to the name
#' of the focal species (set to `NULL` to ignore).
#' @param layout Character string specifying the tree layout passed to
#' \code{ggtree} (e.g. "circular", "rectangular"). Default: `"circular"`.
#' @param offset Number indicating the offset of tip labels from the
#' phylogenetic tree. Default: `5`.
#' @param size Number indicating the size of tip labels. Default: `1.7`.
#' @param tip_colours Character vector of length 2, indicating colours for
#' non-hosts and hosts. Default: `c("black", "coral")`.
#' @param legend_pos Character string specifying legend position. Default:
#' `"none"`.
#'
#' @return A \code{ggplot} object produced by \code{ggtree}.
#'
#' @seealso \code{\link{get_host_spp}}
#'
#' @examples
#' \dontrun{
#' plot_phylogeny(phylo = oak_phylo,
#'                interaction_df = interaction_data,
#'                species = "Agrilus coxalis",
#'                species_col = "agrilus_sp",
#'                host_status_col = "interaction",
#'                host_col = "quercus_sp")
#' }
#'
#' @importFrom ggplot2 aes scale_color_manual theme
#' @importFrom ggtree ggtree geom_tiplab fortify
#' @importFrom dplyr left_join
#' @importFrom treeio as.treedata
#' @export
plot_phylogeny <- function(phylo,
                           host_spp = NULL,
                           interaction_df = NULL,
                           species = NULL,
                           species_col = NULL,
                           host_col = NULL,
                           host_status_col = NULL,
                           main = species,
                           layout = "circular",
                           offset = 5,
                           size = 1.7,
                           tip_colours = c("black", "coral"),
                           legend_pos = "none") {
  # Checks
  if (class(phylo) != "phylo") {
    stop("'phylo' must be an object of class 'phylo'")
  }
  if (!is.null(host_spp) && !is.character(host_spp)) {
    stop("'host_spp' must be a character vector or NULL")
  }
  if (!is.null(!is.character(main)) && !is.character(main)
      || length(main) != 1) {
    stop("'main' must be a single character string or NULL")
  }
  if (!is.character(layout) || length(layout) != 1) {
    stop("'layout' must be a character string of length 1")
  }
  if (!is.numeric(offset) || length(offset) != 1) {
    stop("'offset' must be a single numeric value")
  }
  if (!is.numeric(size) || length(size) != 1) {
    stop("'size' must be a single numeric value")
  }
  if (!is.character(tip_colours) || length(tip_colours) != 2) {
    stop("'tip_colours' must be a character vector of length 2")
  }
  if (!is.character(legend_pos) || length(legend_pos) != 1) {
    stop("'legend_pos' must be a single character string")
  }

  # Get hosts of given species
  if (is.null(host_spp)) {
    host_spp <- get_host_spp(df = interaction_df,
                             species_col = species_col,
                             host_col = host_col,
                             host_status_col = host_status_col,
                             species = species)
  }

  hosts <- setNames(
    data.frame(phylo$tip.label,
               phylo$tip.label %in% host_spp),
    c("taxon", "is_host")
  )

  # Add host info to phylo
  phylo_hosts <- phylo |>
    ggtree::fortify() |>
    dplyr::left_join(hosts, by = c("label" = "taxon")) |>
    treeio::as.treedata()

  # Plot
  ggtree::ggtree(phylo_hosts, layout = layout) +
    ggtree::geom_tiplab(ggplot2::aes(color = is_host),
                        offset = offset,
                        size = size,
                        fontface = "italic") +
    ggplot2::scale_color_manual(values = tip_colours) +
    ggplot2::theme(legend.position = legend_pos) +
    ggplot2::ggtitle(main)
}


#' Plot phylogenetic distances for host species
#'
#' This function creates a barplot showing phylogenetic distances
#' associated with host species for a given focal species.
#'
#' @param interaction_df A data frame containing species–host interaction data.
#' @param species Character string giving the focal species of interest.
#' @param species_col Character string specifying the column containing species
#' names in `interaction_df`.
#' @param host_col Character string specifying the column containing host
#' species names in `interaction_df`.
#' @param host_status_col Character string specifying the column indicating host
#' status in `interaction_df`.
#' @param dist_col Character string specifying the column containing
#' phylogenetic distances in`interaction_df`.
#' @param host_spp A character string specifying the host species. This can be
#' provided as an alternative to inferring the host species via `get_host_spp`
#' when an `interaction_df` is supplied.  Default: `NULL` (not required if
#' `interaction_df` is provided).
#' @param colour Character string specifying the bar colour. Default: `"coral"`.
#' @param cex_names Numeric value controlling the size of axis labels. Default:
#' `0.5`.
#' @param ylab Character string specifying the y-axis label. Default:
#' `"Distance to another host species"`.
#' @param main Character string indicating plot title. Defaults to the name
#' of the focal species (set to `NULL` to ignore).
#'
#' @return Invisibly returns the heights used in the barplot.
#'
#' @seealso \code{\link{get_host_spp}}
#'
#' @examples
#' \dontrun{
#' barplot_host_distance(interaction_df = interaction_data,
#'                       species = "Agrilus coxalis",
#'                       species_col = "agrilus_sp",
#'                       host_col = "quercus_sp",
#'                       host_status_col = "interaction",
#'                       dist_col = "phylo_dist_min")
#' }
#'
#' @export
barplot_host_distance <- function(interaction_df,
                                  species,
                                  species_col,
                                  host_col,
                                  host_status_col,
                                  dist_col,
                                  host_spp = NULL,
                                  colour = "coral",
                                  cex_names = 0.5,
                                  ylab = "Distance to another host species",
                                  main = species) {
  # Checks
  if (!is.data.frame(interaction_df)) {
    stop(paste0("'interaction_df' must be a data.frame"))
  }
  if (!is.character(dist_col) || length(dist_col) != 1) {
    stop("'dist_col' must be a single character string")
  }
  if (!(dist_col %in% colnames(interaction_df))) {
    stop("dist_col must be a column in 'interaction_df'")
  }
  if (!is.numeric(interaction_df[[dist_col]])) {
    stop("'interaction_df[[dist_col]]' must be a numeric vector")
  }
  if (!is.null(host_spp) && !is.character(host_spp)) {
    stop("'host_spp' must be a character vector or NULL")
  }
  if (!is.character(colour)) {
    stop("'colour' must be a character vector")
  }
  if (!is.numeric(cex_names) || length(cex_names) != 1) {
    stop("'cex_names' must be a single number")
  }
  if (!is.character(ylab) || length(ylab) != 1) {
    stop("'ylab' must be a single character string")
  }
  if (!is.null(!is.character(main)) && !is.character(main)
      || length(main) != 1) {
    stop("'main' must be a single character string or NULL")
  }

  # Get hosts of given species
  if (is.null(host_spp)) {
    host_spp <- get_host_spp(df = interaction_df,
                             species_col = species_col,
                             host_col = host_col,
                             host_status_col = host_status_col,
                             species = species)
  }

  # Get distances
  hosts_dist <- subset(interaction_df,
                       interaction_df[[species_col]] == species
                       & interaction_df[[host_col]] %in% host_spp)
  hosts_phylo_dist <- hosts_dist[, c(host_col, dist_col)]

  # Plot
  barplot(height = hosts_phylo_dist[[dist_col]],
          names.arg = hosts_phylo_dist[[host_col]],
          ylab = ylab,
          las = 2,
          cex.names = cex_names,
          col = colour,
          main = species)
}


#' Plot random effect estimates for host species
#'
#' This function extracts and plots random effect estimates from a
#' \code{brmsfit} model for host species associated with a given focal species.
#'
#' @param model A fitted Bayesian model of class \code{brmsfit}.
#' @param host_spp A character string specifying the host species. This can be
#' provided as an alternative to inferring the host species via `get_host_spp`
#' when an `interaction_df` is supplied.  Default: `NULL` (not required if
#' `interaction_df` is provided).
#' @param interaction_df A data frame containing species-host interaction data.
#' Default: `NULL` (not required if `host_spp` is provided).
#' @param species_col A character string indicating the column name in
#' `interaction_df` that contains species names.  Default: `NULL` (not required
#' if `interaction_df` is provided).
#' @param host_col A character string indicating the column name in both `model`
#' and `interaction_df` that contains host species names.
#' @param host_status_col A character string specifying the column in
#' `interaction_df` that indicates the host status.  Default: `NULL` (not
#' required if `interaction_df` is provided).
#' @param species Character string giving the focal species of interest.
#' @param colour Character string specifying the bar colour. Default: `coral`.
#' @param cex_names Numeric value controlling the size of axis labels. Default:
#' `0.5`.
#' @param ylab Character string specifying the y-axis label. Default:
#' `"Random effect estimates"`.
#' @param main Character string indicating plot title. Defaults to the name
#' of the focal species (set to `NULL` to ignore).
#'
#' @return Invisibly returns the random effect estimates used in the barplot.
#'
#' @seealso \code{\link{get_host_spp}}, \code{brms::ranef}
#'
#' @examples
#' \dontrun{
#' plot_ranef(interaction_df = interaction_data,
#'            model = oak_model,
#'            species = "Agrilus coxalis",
#'            species_col = "agrilus_sp",
#'            host_col = "quercus_sp",
#'            host_status_col = "interaction")
#' }
#'
#' @importFrom brms ranef
#' @export
barplot_ranef <- function(model,
                          host_spp = NULL,
                          interaction_df = NULL,
                          species_col = NULL,
                          host_col,
                          host_status_col = NULL,
                          species,
                          colour = "coral",
                          cex_names = 0.5,
                          ylab = "Random effect estimates",
                          main = species) {
  # Checks
  if (!inherits(model, "brmsfit")) {
    stop("'model' must be an object of class brmsfit")
  }
  if (!is.null(host_spp) && !is.character(host_spp)) {
    stop("'host_spp' must be a character vector or NULL")
  }
  if (!is.character(colour)) {
    stop("'colour' must be a character vector")
  }
  if (!is.character(ylab) || length(ylab) != 1) {
    stop("'ylab' must be a single character string")
  }
  if (!is.null(!is.character(main)) && !is.character(main)
      || length(main) != 1) {
    stop("'main' must be a single character string or NULL")
  }

  # Get hosts of given species
  if (is.null(host_spp)) {
    host_spp <- get_host_spp(df = interaction_df,
                             species_col = species_col,
                             host_col = host_col,
                             host_status_col = host_status_col,
                             species = species)
  }

  # Get random effect estimates
  ranef_list <- lapply(host_spp,
                       function(sp) brms::ranef(model)[[host_col]][, , 1][sp, ])
  ranef_df <- do.call(rbind, ranef_list)

  # Plot
  barplot(height = ranef_df[, "Estimate"],
          names.arg = host_spp,
          ylab = ylab,
          las = 2,
          cex.names = cex_names,
          col = colour,
          main = species)
}


#' Extract host species for a given focal species
#'
#' This function identifies host species associated with a focal species
#' from a species–host interaction data frame.
#'
#' @param df A data frame containing species–host interaction data.
#' @param species_col Character string specifying the column containing species
#' names.
#' @param host_col Character string specifying the column containing host
#' species names.
#' @param host_status_col Character string specifying the column indicating
#' host status (values must include 0 and 1, where 1 indicates a host).
#' @param species Character string giving the focal species of interest.
#' @param main Character string indicating plot title. Defaults to the name
#' of the focal species (set to `NULL` to ignore).
#'
#' @return A character vector of unique host species associated with the focal
#' species.
#'
#' @keywords internal
get_host_spp <- function(df,
                         species_col,
                         host_col,
                         host_status_col,
                         species) {
  # Checks
  if (!is.character(species_col) || length(species_col) != 1) {
    stop("'species_col' must be a single character string")
  }
  if (!is.character(host_col) || length(host_col) != 1) {
    stop("'host_col' must be a single character string")
  }
  if (!is.character(host_status_col) || length(host_status_col) != 1) {
    stop("'host_status_col' must be a single character string")
  }

  if (!is.data.frame(df)
      || !(all(c(species_col, host_col, host_status_col) %in% colnames(df)))) {
    stop(paste("'df' must be a data.frame with colnames:",
               species_col, ",", host_status_col, ",", host_col))
  }
  if (!is.character(df[[species_col]])) {
    stop("'df[[species_col]]' must be a character vector")
  }
  if (!is.character(df[[host_col]])) {
    stop("'df[[host_col]]' must be a character vector")
  }
  if (!any(df[[host_status_col]] %in% c(0, 1))) {
    stop("'df[[host_status_col]]' must be a vector of 0s and 1s")
  }
  if (!is.character(species) || length(species) != 1) {
    stop("'species' must be a sigle character string")
  }
  if (!(species %in% df[[species_col]])) {
    stop("'species' must be present in 'df[[species_col]]'")
  }

  # Extract host species
  hosts_df <- unique(subset(df,
                            df[[species_col]] == species
                            & df[[host_status_col]] == 1))
  hosts <- hosts_df[[host_col]]

  hosts
}


#' Plot comparison of cross-validated predictions
#'
#' Creates a scatter plot comparing predictions from two cross-validation
#' methods, with points categorised according to reported status and a
#' classification threshold applied to model predictions.
#'
#' Each observation is classified into one of four groups:
#' \itemize{
#'   \item \code{true_neg}: not reported and predicted below threshold
#'   \item \code{false_pos}: not reported but predicted above threshold
#'   \item \code{false_neg}: reported but predicted below threshold
#'   \item \code{true_pos}: reported and predicted above threshold
#' }
#'
#' These classes are visualised using different point shapes and fills. A
#' reference 1:1 line and horizontal/vertical threshold lines are added to aid
#' interpretation.
#'
#' In addition, Pearson correlation coefficients are computed and displayed on
#' the plot:
#' \itemize{
#'   \item Overall correlation between \code{cv_method_1} and \code{cv_method_2}
#'   \item Correlation restricted to observations with reported positive status
#' }
#'
#' @param cv_method_1 Numeric vector. Predictions from the first
#' cross-validation method.
#' @param cv_method_2 Numeric vector. Predictions from the second
#' cross-validation method. Must be the same length as \code{cv_method_1}.
#' @param reported_status Numeric, factor, or logical vector indicating observed
#' status for each observation. Values equal to 1 are treated as "reported
#' positive", values equal to 0 as "reported negative".
#' @param model_predictions Numeric vector of model predictions used for
#' threshold-based classification. Must be the same length as
#' \code{cv_method_1}.
#' @param threshold Numeric. Threshold on the log-odds scale used to classify
#' predictions as positive or negative.
#' @param xlab Character string. Label for the x-axis.
#' @param ylab Character string. Label for the y-axis.
#'
#' @return A \code{ggplot2} plot object.
#'
#' @details
#' Classification is performed using \code{model_predictions} and
#' \code{threshold}.
#'
#' @examples
#' \dontrun{
#' plot_cv_comparison(
#'   cv_method_1 = loo_lodds,
#'   cv_method_2 = kfold_lodds,
#'   reported_status = host_status,
#'   model_predictions = prediction_lodds,
#'   threshold = threshold
#' )
#' }
#'
#' @import ggplot2
#' @export
plot_cv_comparison <- function(cv_method_1,
                               cv_method_2,
                               reported_status,
                               model_predictions,
                               threshold,
                               xlab = "CV method 1 predictions",
                               ylab = "CV method 2 predictions") {

  # Checks
  if (!is.numeric(cv_method_1)) {
    stop("'cv_method_1' must be a numeric vector")
  }
  if (!is.numeric(cv_method_2)) {
    stop("'cv_method_2' must be a numeric vector")
  }
  if (length(cv_method_2) != length(cv_method_1)) {
    stop("'cv_method_1' and 'cv_method_2' must have the same length")
  }
  if (!(is.logical(reported_status)
        || is.factor(reported_status)
        || is.numeric(reported_status)
        || all(reported_status %in% c(0, 1)))) {
    stop("'reported_status' must be logical, factor, or binary (0/1) numeric")
  }
  if (length(reported_status) != length(cv_method_1)) {
    stop("'reported_status' must match length of 'cv_method_1")
  }
  if (!is.numeric(model_predictions)) {
    stop("'model_predictions' must be a numeric vector")
  }
  if (length(model_predictions) != length(cv_method_1)) {
    stop("'model_predictions' must match length of 'cv_method_1")
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("'model_predictions' must be a single numeric value")
  }
  if (!is.character(xlab) || length(xlab) != 1) {
    stop("'xlab' must be single character string")
  }
  if (!is.character(ylab) || length(ylab) != 1) {
    stop("'ylab' must be single character string")
  }

  # Build plotting data frame
  df <- data.frame(
    cv_method_1 = cv_method_1,
    cv_method_2 = cv_method_2,
    reported_status = reported_status,
    model_predictions = model_predictions
  )

  # Classify interactions
  df$class <- with(df,
                   ifelse(reported_status != 1
                          & model_predictions < threshold,
                          "true_neg",
                          ifelse(reported_status != 1
                                 & model_predictions > threshold,
                                 "false_pos",
                                 ifelse(reported_status == 1
                                        & model_predictions < threshold,
                                        "false_neg",
                                        "true_pos"))))

  df$class <- factor(df$class,
                     levels = c("true_neg", "false_pos",
                                "false_neg", "true_pos"))

  # Compute correlations
  cor_all <- cor(cv_method_1, cv_method_2)
  cor_pos <- cor(cv_method_1[which(reported_status == 1)],
                 cv_method_2[which(reported_status == 1)])

  # Plot
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = cv_method_1,
                                    y = cv_method_2,
                                    shape = class,
                                    fill = class,
                                    alpha = class)) +

    ggplot2::geom_point(colour = "black", size = 2) +

    ggplot2::scale_shape_manual(values = c(
      true_neg = 16,   # filled circle
      false_pos = 17,  # filled triangle
      false_neg = 21,  # filled circle (border)
      true_pos = 24    # filled triangle (border)
    )) +

    ggplot2::scale_fill_manual(values = c(
      true_neg = "black",
      false_pos = "black",
      false_neg = "steelblue",
      true_pos = "steelblue"
    )) +

    ggplot2::scale_alpha_manual(values = c(
      true_neg = 0.05,
      false_pos = 0.05,
      false_neg = 1,
      true_pos = 1
    )) +

    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed") +

    ggplot2::annotate(
      "text",
      x = min(df$cv_method_1),
      y = max(df$cv_method_2),
      hjust = 0,
      vjust = 1,
      label = paste0(
        "Pearson's r (all) = ", round(cor_all, 3), "\n",
        "Pearson's r (known positives) = ", round(cor_pos, 3)
      )
    ) +

    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  p
}
