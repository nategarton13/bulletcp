## end to end data preprocessing and gibbs algorithm
detect_cp_v2 <- function(data, iter = 10000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol_edge = 50, tol_cp = 1000, warmup = 5000, verbose = FALSE,
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
  if(which.max(test_variable_cp_gibbs$max_lp) == 1)
  {
    grooves <- range(nud$x)
  }
  if(which.max(test_variable_cp_gibbs$max_lp) %in% c(2,3))
  {
    upper <- max(d$x[!is.na(d$y)])
    lower <- min(d$x[!is.na(d$y)])
    tmp <- c(test_variable_cp_gibbs$cp_map[[2]]$left, upper)
    tmp2 <- c(lower, test_variable_cp_gibbs$cp_map[[2]]$right)
    if(test_variable_cp_gibbs$max_lp$cp1_left > test_variable_cp_gibbs$max_lp$cp1_right)
    {
      grooves <- tmp
    }
    else{grooves <- tmp2}
  }
  if(which.max(test_variable_cp_gibbs$max_lp) == 4)
  {
    upper <- max(d$x[!is.na(d$y)])
    lower <- min(d$x[!is.na(d$y)])
    tmp <- c(max(lower, test_variable_cp_gibbs$cp_map[[1]][1]), min(upper, test_variable_cp_gibbs$cp_map[[1]][2]))
    grooves <- tmp
  }
  return(list("changepoint_results" = test_variable_cp_gibbs, "cutoffs" = range(nud$x), "grooves" = grooves))
}

