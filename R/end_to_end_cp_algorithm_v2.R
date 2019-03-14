#'  Impute data and estimate groove locations.
#'
#' This function is mostly just a wrapper function which calls the functions necessary
#' to impute missing data, run the changepoint Gibbs algorithms, and select MAP estimates
#' of the changepoint locations. Much less output is given for this function than for the
#' functions called by this function. If all goes well, one should only need to explicitly
#' use this function to estimate groove locations. Note that because this function calls the
#' functions which do the Gibbs sampling, all of the input required for those functions is
#' required by this function.
#' @param data Data frame with columns "x" and "y." "x" is a column of the locations of the
#' observed residual values, y.
#' @param iter Number of iterations after warmup.
#' @param start.vals Starting values for the changepoint algorithm. Either NA valued or a named list
#' of lists. If list, the names of the lists should be "cp2","cp1", and "cp0". Each list posessing
#' one of those aforementioned names is a list of starting values identical to what would be given
#' if the changepoint algorithm were to be run with the corresponding number of specified changepoints.
#' List with elements "sigma", "l", "cp", "beta", and "intercept." "sigma" and "l"
#'   are 3 element vectors where the first element is for the data on the left groove.
#'   The second element is for the land engraved area, and the third element is for the right groove.
#'   "cp" is the vector of changepoint starting values. "beta" and "intercept" are two element vectors
#'   of the slope and intercept for the left and right groove engraved area respectively. If NA,
#'   default starting values will be used. Note that the changepoint starting values should always be
#'   near the edges of the data.
#' @param prop_var Either NA valued or a list of named lists. If list, the names of the lists should be "cp2","cp1", and "cp0".
#' Each list posessing
#' one of those aforementioned names is a list of proposal covariance matrices identical to what would be given
#' if the changepoint algorithm were to be run with the corresponding number of specified changepoints.
#' @param cp_prop_var The proposal variance-covariance matrix for the changepoints. Can either be
#' NA or a named list. If list, the names of the list items should be "cp2", "cp1" where each is the appropriate
#' proposal variance/covariance matrix for the number of changepoints.
#' @param tol_edge This parameter controls how close changepoint proposals can be to the edge of the data
#' before getting automatically rejected. For example, a value of 10 means that the changepoint will be
#' automatically rejected if either of the proposal changepoints is within a distance of 10 x-values from either edge.
#' @param tol_cp This parameter controls how close changepoint proposals can be to each other
#' before getting automatically rejected. For example, a value of 10 means that the changepoint will be
#' automatically rejected if either of the proposal changepoints is within a distance of 10 x-values from either each other.
#' @param warmup The number of warmup iterations. This should be set to a very small number of iterations,
#' as using too many iterations as warmup risks moving past the changepoints and getting stuck in a local mode.
#' Default is set to 500.
#' @param verbose Logical value indicating whether to print the iteration number and the parameter proposals.
#' @param prior_numcp This is a vector with four elements giving the prior probabilities for the zero changepoint model,
#' the one changepoint on the left model, the one changepoint on the right model, and the two changepoint model, in that order.
#' Note that, practically, because the likelihood values are so large, only very strong priors will influence the results.
#' @param est_impute_par Logical value indicating whether parameters for the Gaussian process
#' imputation should be estimated before actually doing the imputation. Default is FALSE, in which case
#' the default imputation standard deviation is 0.8 and the length scale is 15. The covariance function
#' is a squared exponential. These values have worked well in testing.
#' @param impute_par A two element vector containing the standard deviation and length scale (in that order)
#' to use for the Gaussian process imputation. These values will not be used if the est_impute_par
#' argument is set to TRUE.
#' @return A named list containing the output from variable_cp_gibbs function, the range of
#' data that was actually used for the changepoint algorithm (since it doesn't impute values
#' past the outermost non-missing values), and the estimated groove locations.
#' @importFrom stats complete.cases
#' @importFrom Rdpack reprompt
#' @examples
#' # Fake data
#' sim_groove <- function(beta = c(-0.28,0.28), a = 125)
#' {
#'     x <- seq(from = 0, to = 2158, by = 20)
#'     med <- median(x)
#'     y <- 1*(x <= a)*(beta[1]*(x - med) - beta[1]*(a - med)) +
#'     1*(x >= 2158 - a)*(beta[2]*(x - med) - beta[2]*(2158 - a - med))
#'     return(data.frame("x" = x, "y" = y))
#' }
#'
#' fake_groove <- sim_groove()
#' cp_gibbs2 <- detect_cp(data = fake_groove,
#'                     verbose = FALSE,
#'                     tol_edge = 50,
#'                     tol_cp = 1000,
#'                     iter = 300,
#'                     warmup = 100,
#'                     est_impute_par = FALSE)
#'
#' @export

detect_cp <- function(data, iter = 5000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol_edge = 50, tol_cp = 1000, warmup = 200, verbose = FALSE,
                         prior_numcp = rep(1/4, times = 4), est_impute_par = FALSE, impute_par = c(0.8,15))
{
  ## put extra functions in here just in case
  lognormal_ou_pdf <- function(x, mu, sigma, l)
  {
    n <- length(x)
    rho <- exp(-1/l)

    return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
           - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                                - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
  }

  ######### end extra functions

  ## the line immediately below should not be necessary
  # ## this is because hamby44_eval$ccdata_w_resid is actually a list with one element
  # data <- data[[1]]

  ## make range of x data equal to the range of non NA points
  d <- data.frame("x" = data$x, "y" = scale(data$y))
  max_x_notNA <- max(d$x[!is.na(d$y)])
  min_x_notNA <- min(d$x[!is.na(d$y)])
  d <- d[d$x >= min_x_notNA & d$x <= max_x_notNA,]

  ## Put missing values in if there are irregularly spaced gaps in x
  x_na <- seq(from = min(d$x), to = max(d$x), by = min(d$x[2:nrow(d)] - d$x[1:(nrow(d) - 1)]))
  x_na <- x_na[!round(x_na, digits = 2) %in% round(d$x, digits = 2)]
  y_na <- rep(NA, times = length(x_na))
  d_na <- data.frame("x" = x_na, "y" = y_na)
  d <- rbind(d, d_na)
  d <- d[order(d$x),]

  ## put in check to make sure that the tol_edge argument is large enough
  while(nrow(d[d$x > (max(d$x) - tol_edge),]) < 2 ||
        nrow(d[d$x < (min(d$x) + tol_edge),]) < 2)
  {
    tol_edge <- tol_edge + 5
  }

  ## impute data
  temp_d <- d[seq(from = 1, to = nrow(d),by = 20),]
  if(est_impute_par == TRUE)
  {
    temp_dnarm <- temp_d[complete.cases(temp_d),]
    mles <- mlgp(y = temp_dnarm$y, x = temp_dnarm$x)
    impute_par <- exp(mles$par)
  }
  if(any(is.na(d$y)) == TRUE)
  {
    nud <- imputeGP(y = d$y, x = d$x, sigma = impute_par[1], l = impute_par[2])
  }
  else{nud <- d}

  ## run cp algorithm
  test_variable_cp_gibbs <- runmcmc_cpall(data = nud,
                                                          start.vals = start.vals,
                                                          prop_var = prop_var, cp_prop_var = cp_prop_var, verbose = FALSE, tol_edge = tol_edge, tol_cp = tol_cp, iter = iter, warmup = warmup, prior_numcp = prior_numcp)
  ## note that the order that the MAPs are returned in
  ## is 2cp, 1cp left , 1cp right
  if(which.max(test_variable_cp_gibbs$max_lpost) == 1)
  {
    grooves <- range(nud$x)
  }
  if(which.max(test_variable_cp_gibbs$max_lpost) %in% c(2,3))
  {
    upper <- max(d$x[!is.na(d$y)])
    lower <- min(d$x[!is.na(d$y)])
    tmp <- c(test_variable_cp_gibbs$cp_map[[2]]$left, upper)
    tmp2 <- c(lower, test_variable_cp_gibbs$cp_map[[2]]$right)

    if(test_variable_cp_gibbs$max_lpost$cp1_left > test_variable_cp_gibbs$max_lpost$cp1_right)
    {
      grooves <- tmp
    }
    else{grooves <- tmp2}
  }
  if(which.max(test_variable_cp_gibbs$max_lpost) == 4)
  {
    upper <- max(d$x[!is.na(d$y)])
    lower <- min(d$x[!is.na(d$y)])
    tmp <- c(max(lower, test_variable_cp_gibbs$cp_map[[1]][1]), min(upper, test_variable_cp_gibbs$cp_map[[1]][2]))
    grooves <- tmp
  }
  return(list("changepoint_results" = test_variable_cp_gibbs, "cutoffs" = range(nud$x), "grooves" = grooves))
}

#'  Conforming get_grooves_"name" function.
#'
#' This is a wrapper function that comforms to the other get_grooves functions.
#'
#' @param x numeric vector of locations in microns
#' @param value numeric vector of surface measurements in microns
#' @param adjust positive number to adjust the grooves - XXX should be
#'          expressed in microns rather than an index
#' @param ... Additional arguments to be passed to detect_cp_v2.
#' @return A named list containing the output from variable_cp_gibbs function, the range of
#' data that was actually used for the changepoint algorithm (since it doesn't impute values
#' past the outermost non-missing values), and the estimated groove locations.
#' @importFrom stats complete.cases
#' @importFrom stats predict
#' @importFrom dplyr mutate
#' @importFrom stats loess
#' @importFrom assertthat assert_that
#' @importFrom assertthat has_name
#' @examples
#' data("example_data")
#' head(raw_data)
#' raw_data <- raw_data[seq(from = 1, to = nrow(raw_data), by = 30),]
#' cp_gibbs3 <- get_grooves_bcp(x = raw_data$x,
#'     value = raw_data$value,
#'     adjust = 10,
#'     iter = 300,
#'     warmup = 100)
#' @export

get_grooves_bcp <- function(x, value, adjust = 10, ...)
{
  ## get robust loess residuals
  land <- data.frame(x = x, value = value)
  original_land <- land

  ## generate additional variables

  check_min <- min(land$value[!is.na(land$value)])
  land <- dplyr::mutate(land, value_std = value - check_min)
  #install.packages("locfit")
  #library(locfit)

  ## Kiegan's/Susan's robust loess fit function
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

  rlo_fit <- robust_loess_fit(cc = land, iter = 20)
  land$rlo_pred <- predict(rlo_fit, newdata = land)
  land$rlo_resid <- land$value_std - land$rlo_pred

  ## old
  # robust_loess_fit <- locfit::locfit.robust(value_std~x, data = land, alpha = 1, kern = "tcub")
  # land$rlo_pred <- stats::predict(robust_loess_fit, newdata = land)
  # land$rlo_resid <- with(land, value_std-rlo_pred)

  ## create data frame to be passed to detect_cp_v2
  data <- data.frame("x" = land$x, "y" = land$rlo_resid)
  # data <- data.frame("x" = land$x, "y" = land$value)

  cp_results <- detect_cp(data = data, ...)

  ## groove locations plus adjustment
  groove <- c(cp_results$grooves[1] + adjust, cp_results$grooves[2] - adjust)


  return(list(groove = groove, "changepoint_samples" = cp_results$changepoint_results))

}
