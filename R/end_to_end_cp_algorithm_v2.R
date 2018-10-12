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
#' @param prior_numcp Vector of three values giving the prior weights for the number of grooves/changepoints.
#' The default value is a uniform prior.
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
#' @export

detect_cp_v2 <- function(data, iter = 5000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol_edge = 50, tol_cp = 1000, warmup = 500, verbose = FALSE,
                         prior_numcp = c(1/3, 1/3, 1/3), est_impute_par = FALSE, impute_par = c(0.8,15))
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

  ## this is because hamby44_eval$ccdata_w_resid is actually a list with one element
  data <- data[[1]]

  ## make range of x data equal to the range of non NA points
  d <- data.frame("x" = data$y, "y" = scale(data$rlo_resid))
  max_x_notNA <- max(d$x[!is.na(d$y)])
  min_x_notNA <- min(d$x[!is.na(d$y)])
  d <- d[d$x >= min_x_notNA & d$x <= max_x_notNA,]

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
    nud <- myimpute(y = d$y, x = d$x, sigma = impute_par[1], l = impute_par[2])
  }
  else{nud <- d}

  ## run cp algorithm
  test_variable_cp_gibbs <- variable_cp_gibbs_v2(data = nud,
                                                          start.vals = start.vals,
                                                          prop_var = prop_var, cp_prop_var = cp_prop_var, verbose = FALSE, tol_edge = tol_edge, tol_cp = tol_cp, iter = iter, warmup = warmup, prior_numcp = prior_numcp)
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

