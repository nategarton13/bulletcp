## end to end data preprocessing and gibbs algorithm
detect_cp_v2 <- function(data, iter = 10000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol = 10, warmup = 5000, verbose = FALSE,
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

  ## remove NA values from the data
  d <- data.frame("x" = data$y, "y" = scale(data$rlo_resid))

  ## impute data
  temp_d <- d[seq(from = 1, to = nrow(d),by = 20),]
  if(est_impute_par == TRUE)
  {
    temp_dnarm <- temp_d[complete.cases(temp_d),]
    mles <- mlgp(y = temp_dnarm$y, x = temp_dnarm$x)
    impute_par <- exp(mles$par)
  }
  nud <- myimpute(y = d$y, x = d$x, sigma = impute_par[1], l = impute_par[2])

  ## run cp algorithm
  test_variable_cp_gibbs <- variable_cp_gibbs_v2(data = nud,
                                                          start.vals = start.vals,
                                                          prop_var = prop_var, cp_prop_var = cp_prop_var, verbose = FALSE, tol = tol, iter = iter, warmup = warmup, prior_numcp = prior_numcp)
  if(which.max(test_variable_cp_gibbs$avg_lp) == 1)
  {
    grooves <- range(nud$x)
  }
  if(which.max(test_variable_cp_gibbs$avg_lp) == 2)
  {
    upper <- max(d$x[!is.na(d$y)])
    lower <- min(d$x[!is.na(d$y)])
    tmp <- c(test_variable_cp_gibbs$cp_mean[[2]]$left, upper)
    tmp2 <- c(lower, test_variable_cp_gibbs$cp_mean[[2]]$right)
    if(test_variable_cp_gibbs$avg_lp_1cp[1] > test_variable_cp_gibbs$avg_lp_1cp[2])
    {
      grooves <- tmp
    }
    else{grooves <- tmp2}
  }
  if(which.max(test_variable_cp_gibbs$avg_lp) == 3)
  {
    grooves <- test_variable_cp_gibbs$cp_mean[[1]]
  }
  return(list("changepoint_results" = test_variable_cp_gibbs, "cutoffs" = range(nud$x), "grooves" = grooves))
}

