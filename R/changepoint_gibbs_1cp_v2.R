#'  Estimate a posterior distribution of data conditional that there is one groove.
#'
#' This function is basically a wrapper for running the left and right (one) changepoint
#' Gibbs algorithms. The only computation that this function does is to estimate the posterior
#' means of the left and right changepoint distributions.
#' @param data Data frame with columns "x" and "y." "x" is a column of the locations of the
#' observed residual values, y.
#' @param iter Number of interations after warmup.
#' @param start.vals.left Starting values for the changepoint algorithm assuming the groove is on the left.
#' List with elements "sigma", "l", "cp", "beta", and "intercept." "sigma" and "l"
#'   are 2 element vectors where the first element is for the data to the left of the changepoint.
#'   "cp" is the changepoint starting value. "beta" and "intercept" are the slope and intercept
#'   starting values for the mean of the data model to the left of the changepoint.
#'  which parameterize the covariance matrix.
#' @param start.vals.right Starting values for the changepoint algorithm assuming the groove is on the right.
#' @param prop_var_left The proposal variance for the random walk Metropolis algorithm assuming
#' that the groove is on the left.
#' A two element list of the proposal variance-covariance matrices for the random
#' walk metropolis algorithm(s). The first element is for the data to the left of the changepoint.
#' @param prop_var_right The proposal variance for the random walk Metropolis algorithm assuming
#' that the groove is on the right.
#' @param cp_prop_var The proposal variance for the changepoint.
#' @param tol_edge This parameter controls how close changepoint proposals can be to the edge of the data
#' before getting automatically rejected. For example, a value of 10 means that the changepoint will be
#' automatically rejected if the proposal is within a distance of 10 x-values from either edge.
#' @param warmup The number of initial iterations which serves two purposes: the first is to allow the
#' algorithm to wander to the area of most mass, and the second is to tune the proposal variance.
#' @param verbose Logical value indicating whether to print the iteration number and the parameter proposals.
#' @return A named list with all of the output that the left and right changepoint functions produce
#' individually plus the posterior means of the left and right changepoints.
#' @export
runmcmc_cp1 <- function(data, iter, start.vals.left, start.vals.right, prop_var_left, prop_var_right, cp_prop_var, tol_edge = 10, warmup = 500, verbose = FALSE)
{
  ##data is a data frame with column x and column y

  ## run left cp algorithm first
  left_cp_out <- runmcmc_cp1_left(data = data, iter = iter, warmup = warmup, start.vals = start.vals.left, prop_var = prop_var_left, cp_prop_var = cp_prop_var, verbose = verbose, tol_edge = tol_edge)

  ## run right cp algorithm
  right_cp_out <- runmcmc_cp1_right(data = data, iter = iter, warmup = warmup, start.vals = start.vals.right, prop_var = prop_var_right, cp_prop_var = cp_prop_var, verbose = verbose, tol_edge = tol_edge)

  ## compute posterior probabilities of left or right changepoint
  # prior <- 0.5
  # pleft <- (prior * mean(left_cp_out$lp)) / (prior * mean(left_cp_out$lp) + (1-prior) * mean(right_cp_out$lp))
  # pright <- 1 - pleft
  # pleft <- NA
  # pright <- NA

  ## estimated changepoint
  est_left_cp <- mean(left_cp_out$parameters$cp)
  est_right_cp <- mean(right_cp_out$parameters$cp)

  return(list("left_parameters" = left_cp_out$parameters, "right_parameters" = right_cp_out$parameters,
              "accept" = list("left" = left_cp_out$accept, "right" = right_cp_out$accept),
              "lp" = list("left" = left_cp_out$lp, "right" = right_cp_out$lp),
              "lpost" = list("left" = left_cp_out$lpost, "right" = right_cp_out$lpost),
              "estimated_cp" = list("left" = est_left_cp, "right" = est_right_cp))
         )
}

