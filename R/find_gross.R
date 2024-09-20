#' @title Find the gross outliers.
#'
#' @inheritParams find_elbow
#' @inheritParams outcast_gmm
#' @param k_neighbours Number of neighbours for dbscan::glosh
#' @param choice Optional preset number of gross outliers.
#'
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $bool: a logical vector identifying the gross outliers.
#' * $plot
#'
#' @export
find_gross <- function(
    x, search_centre,
    k_neighbours = floor(nrow(x) / 100),
    underestimate = 0.5,
    choice = NULL) {
  
  stopifnot(!is.null(search_centre))
  
  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  knndist_sort <- -sort(-x_knndist)
  
  elbow <- find_elbow(knndist_sort, search_centre, TRUE)

  if (is.null(choice)) {
    choice <- elbow$choice
  }
  
  gross_choice <- floor(choice * underestimate)

  bool <- rank(-x_knndist) <= gross_choice
  
  plot <- elbow$plot + 
    ggplot2::geom_vline(xintercept = gross_choice, colour = "red") +
    ggplot2::labs(
      subtitle = paste0(
        "Gross choice = ", gross_choice
      ),
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    )

  return(list(choice = gross_choice, bool = bool, plot = plot))
}

#' Plot sorted dbscan::glosh values to determine the `search_centre` for
#' `find_gross`.
#'
#' @inheritParams find_gross
#'
#' @return A ggplot2 object
#'
#' @export
plot_gross <- function(x, k_neighbours = floor(nrow(x) / 100)) {

  x_knndist <- dbscan::kNNdist(x, k_neighbours)
  
  outlier_seq <- seq_len(nrow(x))
  knndist_sort <- -sort(-x_knndist)
  gg <- data.frame(outlier_seq, knndist_sort) |>
    ggplot2::ggplot(ggplot2::aes(x = outlier_seq, y = knndist_sort)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Outlier Number",
      y = paste0("kNN Distance (k = ", k_neighbours, ")")
    )

  return(gg)
}

#' @title Find the elbow using a two-piece linear regression model.
#'
#' @param y Numeric vector.
#' @param search_centre Centre of changepoint search interval.
#' @param concave Logical value.
#'
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $search_interval: a numeric vector of length 2 indicating the minimum and
#'                     maximum values considered for the changepoint.
#'
#' @export
find_elbow <- function(y, search_centre = NULL, concave = TRUE) {
  y_len <- length(y)
  outliers_removed <- seq_len(y_len)
  
  if (is.null(search_centre)) {
    linmod <- stats::lm(y ~ seq_len(y_len))
    if (concave) {
      search_centre <- which.min(linmod$residuals)
    } else {
      search_centre <- which.max(linmod$residuals)
    }
  }
  
  search_radius <- min(c(
    floor((y_len - search_centre - 1) / 3),
    floor((search_centre - 2) / 3)
  ))
  
  if (search_radius == 0 && search_centre < 5) {
    warning("search_radius == 0. Try search_centre > 4.\n")
  } else if (search_radius == 0 && search_centre > (y_len - 4)) {
    warning(paste0(
      "search_radius == 0. Try search_centre < ", y_len - 3, ".\n"
    ))
  }
  
  upper <- search_centre + search_radius
  lower <- search_centre - search_radius
  search_interval <- c(lower, upper)
  
  width <- floor(y_len / 10)
  left <- right <- integer(y_len)
  sd_vec <- double(y_len)
  for (i in outliers_removed) {
    left[i] <- min(max(c(i - width, 1)), y_len - 2 * width)
    right[i] <- max(min(c(i + width, y_len)), 1 + 2 * width)
    sd_vec[i] <- sqrt(mean(diff(diff(y[left[i]:right[i]]))^2) / 6)
  }
  
  cpop_out <- suppressMessages(cpop::cpop(
    y, outliers_removed,
    grid = lower:upper, minseglen = 2 * search_radius + 1,
    sd = sd_vec
  ))
  
  stopifnot(length(cpop_out@changepoints) == 3)
  
  choice <- cpop_out@changepoints[2]
  
  cat("search centre = ", search_centre)
  cat(", search radius = ", search_radius)
  cat(", search interval = ", search_interval)
  cat(", search choice = ", choice, "\n")
  
  cpop_fitted <- cpop::fitted(cpop_out)
  
  gg <- cbind(x = outliers_removed, y = y) |>
    ggplot2::ggplot(ggplot2::aes(x = outliers_removed, y = y)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = search_interval, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = choice) +
    ggplot2::geom_abline(
      slope = cpop_fitted$gradient[1], intercept = cpop_fitted$intercept[1],
      linetype = "dotted"
    ) +
    ggplot2::geom_abline(
      slope = cpop_fitted$gradient[2], intercept = cpop_fitted$intercept[2],
      linetype = "dotted"
    ) +
    ggplot2::labs(
      title = paste0(
        "Elbow choice = ", choice,
        " (Search interval = [", search_interval[1],
        ", ", search_interval[2], "])"
      ),
      x = "Outlier Number"
    )
  
  return(list(
    choice = choice,
    search_interval = search_interval,
    cpop_out = cpop_out,
    plot = gg
  ))
}
