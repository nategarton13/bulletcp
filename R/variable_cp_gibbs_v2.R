#'  Estimate posterior distributions for the 0, 1, or 2 changepoint case.
#'
#' This function runs the changepoint functions designed for the cases when there are
#' 0, 1, or 2 changepoints. It then returns a subset of the results that are returned for
#' each function individually. This subset of results is enough to decide the likely number of
#' shoulders, the locations of the shoulders (if they exist), as well as the posterior
#' samples for the changepoints for minimal diagnostic use.
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
#' @param prior_numcp This is a vector with four elements giving the prior probabilities for the zero changepoint model,
#' the one changepoint on the left model, the one changepoint on the right model, and the two changepoint model, in that order.
#' Note that, practically, because the likelihood values are so large, only very strong priors will influence the results.
#' @param warmup The number of warmup iterations. This should be set to a very small number of iterations,
#' as using too many iterations as warmup risks moving past the changepoints and getting stuck in a local mode.
#' Default is set to 500.
#' @param verbose Logical value indicating whether to print the iteration number and the parameter proposals.
#' @return A named list containing the sampled changepoint locations for both the one and two changepoint scenarios,
#' the posterior changepoint means, the average log pdf values from the data model under each model,
#' the maximum log probability values under each model
#' log likelihood values, and estimates of the maximum a posteriori changepoint value
#' under each model.
#' @example
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
#'
#' # define starting values for the changepoints
#' cp_start_left <- min(fake_groove$x) + 60
#' cp_start_right <- max(fake_groove$x) - 60
#'
#' # define list of starting values for both the left and right changepoint models
#' cp0.start.vals <- list("sigma" = c(1), "l" = c(10))
#' cp1.start.vals <- list("left" = list("sigma" = c(1,1),
#'                               "l" = c(10,10),
#'                               "cp" = c(cp_start_left),
#'                               "beta" = c(-1),
#'                               "intercept" = c(0)),
#'                               "right" = list("sigma" = c(1,1),
#'                                "l" = c(10,10),
#'                                 "cp" = c(cp_start_right),
#'                                 "beta" = c(1),
#'                                 "intercept" = c(0)))
#' cp2.start.vals <- list("sigma" = c(1,1,1),
#'                 "l" = c(10,10,10),
#'                 "cp" = c(cp_start_left, cp_start_right),
#'                 "beta" = c(-2,2),
#'                 "intercept" = c(0,0))
#' start.vals <- list("cp2" = cp2.start.vals, "cp1" = cp1.start.vals, "cp0" = cp0.start.vals)
#'
#' # list of starting values for each of the two MH steps (not sampling the changepoint) for both the left and right changepoint models
#' prop_var_0cp <- diag(c(1/2,1/2))
#' prop_var_lrcp <- list("left" = list(diag(c(1/2,1/2,1/2,1/2)),
#'                             diag(c(1/2,1/2))),
#'                             "right" = list(diag(c(1/2,1/2)),
#'                             diag(c(1/2,1/2,1/2, 1/2))))
#'
#' prop_var_2cp <- list(diag(c(1/2,1/2,1/2,1/2)),
#'               diag(c(1/2,1/2)),
#'               diag(c(1/2,1/2,1/2,1/2)))
#'
#' prop_var <- list("cp2" = prop_var_2cp, "cp1" = prop_var_lrcp, "cp0" = prop_var_0cp)
#'
#' # define the proposal variance for the RWMH step sampling the changepoint
#' cp_prop_var <- list("cp2" = diag(c(10^2, 10^2)),
#'                  "cp1" = 10^2)
#'
#' # prior on the number of changepoints
#' prior_numcp <- rep(1/4, times = 4)
#'
#' set.seed(1111)
#' cp_gibbs <- runmcmc_cpall(data = fake_groove,
#'                        start.vals = start.vals,
#'                        prop_var = prop_var,
#'                        cp_prop_var = cp_prop_var,
#'                        verbose = FALSE,
#'                        tol_edge = 50,
#'                        tol_cp = 1000,
#'                        iter = 300,
#'                        warmup = 100,
#'                        prior_numcp = prior_numcp)
#' @export


## function to get the conditional posterior given 0,1,2 changepoints
runmcmc_cpall <- function(data, iter = 8000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol_edge = 50, tol_cp = 1000, warmup = 500, verbose = FALSE, prior_numcp = rep(1/4, times = 4))
{
  ## If some function arguments (starting values/proposal variances are unspecified)
  ## choose generic arguments
  if(any(is.na(start.vals)))
  {
    cp_sval_left <- min(data$x) + tol_edge + 10
    cp_sval_right <- max(data$x) - tol_edge - 10
    start.vals <- list("cp2" = list("sigma" = c(1,1,1), "l" = c(10,10,10), "cp" = c(cp_sval_left,cp_sval_right), "beta" = c(-2,2), "intercept" = c(0,0)),
                       "cp1" = list("left" = list("sigma" = c(1,1), "l" = c(10,10), "cp" = c(cp_sval_left), "beta" = c(-1), "intercept" = c(0)),
                                    "right" = list("sigma" = c(1,1), "l" = c(10,10), "cp" = c(cp_sval_right), "beta" = c(1), "intercept" = c(0))),
                       "cp0" = list("sigma" = c(1), "l" = c(10)))
  }
  if(any(is.na(prop_var)))
  {
    prop_var <- list("cp2" = list(diag(c(1/2,1/2,1/2,1/2)), diag(c(1/2,1/2)), diag(c(1/2,1/2,1/2,1/2))),
                     "cp1" = list("left" = list(diag(c(1/2,1/2,1/2,1/2)), diag(c(1/2,1/2))),
                                  "right" = list(diag(c(1/2,1/2)), diag(c(1/2,1/2,1/2, 1/2)))),
                     "cp0" = diag(c(1/2,1/2)))
  }
  if(any(is.na(cp_prop_var)))
  {
    cp_prop_var <- list("cp2" = diag(c(10^2, 10^2)),
                        "cp1" = 10^2)
  }

  ## change point parameter list
  cp_list <- list()

  ## two changepoint model
  cp2_dsn <- runmcmc_cp2(data = data, iter = iter, start.vals = start.vals$cp2, prop_var = prop_var$cp2, cp_prop_var = cp_prop_var$cp2, tol_edge = tol_edge, tol_cp = tol_cp, warmup = warmup, verbose = verbose)
  cp_list$cp2 <- cp2_dsn$parameters$cp
  mcp2 <- mean(exp(cp2_dsn$lp))
  map_cp2 <- cp2_dsn$parameters$cp[which.max(cp2_dsn$lpost),]

  ## one changepoint model
  cp1_dsn <- runmcmc_cp1(data = data, iter = iter, start.vals.left = start.vals$cp1$left, start.vals.right = start.vals$cp1$right,
                          prop_var_left = prop_var$cp1$left, prop_var_right = prop_var$cp1$right, cp_prop_var = cp_prop_var$cp1, tol_edge = tol_edge, warmup = warmup, verbose = verbose)
  ##cp_list$cp1 <- list("ppleft" = cp1_dsn$ppleft, "ppright" = cp1_dsn$ppright, "left" = cp1_dsn$left_parameters$cp, "right" = cp1_dsn$right_parameters$cp)
  cp_list$cp1 <- list("left" = cp1_dsn$left_parameters$cp, "right" = cp1_dsn$right_parameters$cp)
  mcp1 <- 0.5 * mean(exp(cp1_dsn$lp$left)) + 0.5 * mean(exp(cp1_dsn$lp$right))
  map_cp1_left <- as.numeric(cp1_dsn$left_parameters$cp)[which.max(cp1_dsn$lpost$left)]
  map_cp1_right <- as.numeric(cp1_dsn$right_parameters$cp)[which.max(cp1_dsn$lpost$right)]

  ## zero changepoint model
  cp0_dsn <- runmcmc_cp0(data = data, iter = iter, start.vals = start.vals$cp0, prop_var = prop_var$cp0, warmup = warmup, verbose = verbose)
  mcp0 <- mean(cp0_dsn$lp)

  ## posterior cp probabilities
  # ratio02 <- exp(log(prior_numcp[1]) + log(mcp0) - log(prior_numcp[3]) - log(mcp2))
  # ratio12 <- exp(log(prior_numcp[2]) + log(mcp1) - log(prior_numcp[3]) - log(mcp2))
  #
  # p2 <- 1/(ratio02 + ratio12 + 1)
  # p1 <- ratio12 * p2
  # p0 <- ratio02 * p2
  #
  # post_numcp <- c(p0,p1,p2)

  # avg_lp <- list("cp0" = mean(cp0_dsn$lp), "cp1" = 0.5 * mean(cp1_dsn$lp$left) + 0.5 * mean(cp1_dsn$lp$right), "cp2" = mean(cp2_dsn$lp))
  # avg_lp_left <- mean(cp1_dsn$lp$left)
  # avg_lp_right <- mean(cp1_dsn$lp$right)

  max_lp <- list("cp0" = max(cp0_dsn$lp), "cp1_left" = max(cp1_dsn$lp$left), "cp1_right" = max(cp1_dsn$lp$right), "cp2" = max(cp2_dsn$lp))

  max_lpost <- list("cp0" = max(cp0_dsn$lpost) + log(prior_numcp[1]), "cp1_left" = max(cp1_dsn$lpost$left) + log(prior_numcp[2]), "cp1_right" = max(cp1_dsn$lpost$right) + log(prior_numcp[3]), "cp2" = max(cp2_dsn$lpost) + log(prior_numcp[4]))
  cp_map <- list("2cp" = map_cp2, "1cp" = list("left" = map_cp1_left, "right" = map_cp1_right))

  return(list("posterior_cp" = cp_list,
              "cp_mean" = list("2cp" = apply(X = cp_list$cp2, MARGIN = 2, FUN = mean),
                               "1cp" = list("left" = mean(cp_list$cp1$left), "right" = mean(cp_list$cp1$right))),
              "max_lp" = max_lp,
              "max_lpost" = max_lpost,
              "cp_map" = cp_map))

}

# ## test the full procedure
# bullet_resid <- hamby44$ccdata_w_resid[[6]]
# test_dat <- bullet_resid[!is.na(bullet_resid$rlo_resid),]
# d <- data.frame("x" = test_dat$y, "y" = scale(test_dat$rlo_resid))
# plot(d$x, d$y)
#
# ## preprocess data
# temp_d <- preprocess_bullet_resid_simple(data = d)
# points(temp_d$d$x, temp_d$d$y, col = "red")
#
# ## run cp algorithm

lognormal_ou_pdf <- function(x, mu, sigma, l)
{
  n <- length(x)
  rho <- exp(-1/l)

  return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
         - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                              - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
}

## function to sample data weighted according to distance from center
# sample_weights <- function(x, sq = FALSE, offset = 1)
# {
#   center <- (min(x) + max(x))/2
#   ifelse(sq == FALSE, dist <- abs(x - center), dist <- (abs(x - center)^2))
#   weights <- (dist + offset)/sum(dist + offset)
#   return(weights)
# }
