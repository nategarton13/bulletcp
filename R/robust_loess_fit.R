#' Fit a robust loess regression
#'
#' Internal function called by get_grooves_lassobasic and get_grooves_lassofull
#' @param cc data frame with columns x and value_std, representing the crosscut
#' @param iter number of iterations
#' @importFrom stats loess
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
#' @example
#' data("example_data")
#' head(raw_data)
#' raw_data <- raw_data[seq(from = 1, to = nrow(raw_data), by = 30),]
#' plot(raw_data$x, raw_data$y)
#'
#' # set the minimum y-value to zero
#' check_min <- min(raw_data$value[!is.na(raw_data$value)])
#' raw_data <- dplyr::mutate(raw_data, value_std = value - check_min)
#'
#' # remove global structure
#' rlo_fit <- robust_loess_fit(cc = raw_data, iter = 20)
#' raw_data$rlo_pred <- predict(rlo_fit, newdata = raw_data)
#' raw_data$rlo_resid <- raw_data$value_std - raw_data$rlo_pred
#'
#' # define new data frame without the global structure
#' data <- data.frame("x" = raw_data$x, "y" = raw_data$rlo_resid)
#' plot(data$x, data$y)
#' @export
robust_loess_fit <- function(cc, iter) {
  assert_that(has_name(cc, "x"), has_name(cc, "value_std"))
  n <- nrow(cc)
  weights <- rep(1, n)
  fit <- loess(value_std ~ x, data = cc, span = 1)
  cc$fit <- predict(fit, newdata = cc)
  cc$resid <- cc$value_std - cc$fit
  i <- 1
  while (i < iter) {
    mar <- median(abs(cc$resid), na.rm = T)
    cc$bisq <- pmax(1 - (cc$resid / (6 * mar))^2, 0)^2
    weights <- ifelse(cc$resid > 0, cc$bisq, 1)
    fit <- loess(value_std ~ x, data = cc, span = 1, weights = weights)
    cc$fit <- predict(fit, newdata = cc)
    cc$resid <- cc$value_std - cc$fit
    i <- i + 1
  }
  return(fit)
}
