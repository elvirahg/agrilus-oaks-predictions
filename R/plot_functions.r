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
#'   dist_cols = c("geo_min_dist_mean", "geo_min_dist_mean_norm", "geo_min_dist_mean_log")
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
    } else if (colname %in% c("geo_min_dist_mean_norm", "geo_min_dist_median_norm")) {
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
    } else if (colname %in% c("geo_min_dist_mean_norm", "geo_min_dist_median_norm")) {
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
#'                     full = "full_model_pred",
#'                     other = "other_model_pred",
#'                     xlab_text = "Other model prediction")
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
#'   \describe{
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
      stop("'obs_vals' must be a numeric vector of equal length to 'pred_original'")
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

#' Plot phylogeny with predicted and known Hhsts
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
#' pred_hosts <- data.frame(
#'   quercus_sp = oak_phylo$tip.label,
#'   is_known = oak_phylo$tip.label %in% known_hosts,
#'   is_pred = oak_phylo$tip.label %in% pred_hosts
#' )
#' plot_oak_phylo_hosts(oak_phylo, pred_hosts, label = "quercus_sp")
#'
#' @import ggplot2 ggtree treeio
#' @export
plot_oak_phylo_hosts <- function(phylo,
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
      || any(colnames(pred_hosts) != c(label, "is_known", "is_pred"))) {
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
    ggplot2::fortify() |>
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


#' Plot plant species distribution per region
#'
#' This function prepares counts of plant species per region and plots them
#' on a world map with discrete colour bins.
#'
#' @param df A data frame (or sf object) containing plant distribution data.
#' @param group_col Character string specifying the column to group by (default
#' `"LEVEL3_NAM"`).
#' @param counts_type Character string specifying which counts to plot.
#'   Options are `"sp"`, `"native_sp"`, `"hosts"`, `"hosts_native"`, `"pred"`,
#' `"native_pred"`.
#' @param breaks Numeric vector specifying breakpoints for binning counts.
#' @param labels Character vector of labels for the bins.
#' @param world_sf Optional `sf` object with world polygons. If `NULL`, uses
#' rnaturalearth.
#' @param palette_val Colour palette values for `scale_fill_brewer`; must be a
#' character vector. Defaults to RColorBrewer::brewer.pal(9, "PuBuGn").
#' @param na_col Colour for NA values. Defaults to `"gray70"`.
#'
#' @return A ggplot object showing a world map coloured by the selected counts.
#'
#' @examples
#' plot_distribution(oak_distributions, counts_type = "native_sp")
#' @import rnaturalearth ggplot2 RColorBrewer
#' @export
plot_distribution <- function(df,
                              group_col = "LEVEL3_NAM",
                              counts_type,
                              breaks = c(-Inf, 0, 10, 20,
                                         30, 40, Inf),
                              labels = c("0", "1–10", "11–20",
                                         "21–30", "31–40", "40+"),
                              world_sf = NULL,
                              palette_vals = RColorBrewer::brewer.pal(9, "PuBuGn"),
                              na_col = "gray70") {
  # Checks
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.null(world_sf) && !inherits(world, "sf")) {
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
  counts_binned <- prepare_counts(df = df,
                                  group_col = group_col,
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
                                "distribution, grouped by", group_col))

}


#' Prepare counts per region and bin them
#'
#' This is an internal helper function that computes the number of
#' plant species, native species, known hosts, or predicted hosts
#' per region, and assigns each region to a discrete bin for plotting.
#'
#' @param df A data frame (or sf object) containing plant distribution
#' data. Must include the columns `plant_sp`, `occurrence_type`,
#' `known_host` (for host couts), `pred_host` (for predicte host counts),
#' and the grouping column.
#' @param group_col Character string giving the name of the column to group by.
#' Defaults to `"LEVEL3_NAM"`.
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
#'   - the grouping column (`group_col`)
#'   - `value`: the raw count for the selected `counts_type`
#'   - `value_band`: the binned factor for plotting
#'
#' @import dplyr
#' @keywords internal
prepare_counts <- function(df,
                           group_col = "LEVEL3_NAM",
                           counts_type,
                           breaks = c(-Inf, 0, 10, 20, 30, 40, Inf),
                           labels = c("0", "1–10", "11–20", "21–30",
                                      "31–40", "40+")) {
  # Checks
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }
  if (!is.character(group_col) || length(group_col) != 1) {
    stop("'group_col' must be a single character string")
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
    dplyr::group_by(!!dplyr::sym(group_col)) |>
    dplyr::summarise(
      value = dplyr::case_when(
        counts_type == "species" ~ n_distinct(plant_sp),
        counts_type == "species_native" ~ n_distinct(plant_sp[occurrence_type == "native"]),
        counts_type == "hosts" ~ n_distinct(plant_sp[known_host]),
        counts_type == "hosts_native" ~ n_distinct(plant_sp[known_host
                                                            & occurrence_type == "native"]),
        counts_type == "pred_hosts" ~ n_distinct(plant_sp[pred_host]),
        counts_type == "pred_hosts_native" ~ n_distinct(plant_sp[pred_host
                                                                 & occurrence_type == "native"])
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(!!dplyr::sym(group_col)) |>
    dplyr::mutate(value_band = cut(value,
                                   breaks = breaks,
                                   labels = labels))

  cat("Counts value summary:\n")
  summary_result <- summary(counts$value)
  cat(paste(names(summary_result), summary_result, sep = ": "), sep = "\n")

  counts
}
