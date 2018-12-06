#' Fit a robust loess regression
#'
#' Internal function called by get_grooves_lassobasic and get_grooves_lassofull
#' @param cc data frame with columns x and value_std, representing the crosscut
#' @param iter number of iterations
#' @importFrom stats loess
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
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
