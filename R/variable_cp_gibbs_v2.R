lognormal_ou_pdf <- function(x, mu, sigma, l)
{
  n <- length(x)
  rho <- exp(-1/l)

  return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
         - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                              - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
}
lognormal_ou_pdf_par <- function(par, x)
{
  sigma <- par[1]
  l <- par[2]
  mu <- rep(0, times = length(x))
  n <- length(x)
  rho <- exp(-1/l)

  return(-n/2 * log(2 * pi) - n * log(sigma) - ((n - 1)/2) * log(1 - rho^2)
         - 1/2 * 1/(sigma^2 * (1 - rho^2)) * ((x[1] - mu[1])^2 + (x[n] - mu[n])^2 + (1 + rho^2) * sum((x[2:(n-1)] - mu[2:(n-1)])^2)
                                              - 2 * rho * sum((x[1:(n-1)] - mu[1:(n-1)]) * (x[2:n] - mu[2:n]))))
}

## function to sample data weighted according to distance from center
sample_weights <- function(x, sq = FALSE, offset = 1)
{
  center <- (min(x) + max(x))/2
  ifelse(sq == FALSE, dist <- abs(x - center), dist <- (abs(x - center)^2))
  weights <- (dist + offset)/sum(dist + offset)
  return(weights)
}

## function to get the conditional posterior given 0,1,2 changepoints
variable_cp_gibbs_v2 <- function(data, iter = 10000, start.vals = NA, prop_var = NA, cp_prop_var = NA, tol = 10, warmup = 5000, verbose = FALSE, prior_numcp = c(1/3, 1/3, 1/3))
{
  ## If some function arguments (starting values/proposal variances are unspecified)
  ## choose generic arguments
  if(is.na(start.vals))
  {
    start.vals <- list("cp2" = list("sigma" = c(1,1,1), "l" = c(10,10,10), "tau" = c(1,1,1), "cp" = c(1000,1500), "beta" = c(-2,2), "intercept" = c(0,0)),
                       "cp1" = list("left" = list("sigma" = c(1,1), "l" = c(10,10), "tau" = c(1,1), "cp" = c(1000), "beta" = c(-1), "intercept" = c(0)),
                                    "right" = list("sigma" = c(1,1), "l" = c(10,10), "tau" = c(1,1), "cp" = c(1000), "beta" = c(1), "intercept" = c(0))),
                       "cp0" = list("sigma" = c(1), "l" = c(10), "tau" = c(1)))
  }
  if(is.na(prop_var))
  {
    prop_var <- list("cp2" = list(diag(c(1/2,1/2,1/2,1/2,1/2)), diag(c(1/2,1/2,1/2)), diag(c(1/2,1/2,1/2,1/2,1/2))),
                     "cp1" = list("left" = list(diag(c(1/2,1/2,1/2,1/2, 1/2)), diag(c(1/2,1/2,1/2))),
                                  "right" = list(diag(c(1/2,1/2,1/2)), diag(c(1/2,1/2,1/2,1/2, 1/2)))),
                     "cp0" = diag(c(1/2,1/2,1/2)))
  }
  if(is.na(cp_prop_var))
  {
    cp_prop_var <- list("cp2" = diag(c(10^2, 10^2)),
                        "cp1" = 10^2)
  }

  ## change point parameter list
  cp_list <- list()

  ## two changepoint model
  cp2_dsn <- cp2_gibbs_v2(data = data, iter = iter, start.vals = start.vals$cp2, prop_var = prop_var$cp2, cp_prop_var = cp_prop_var$cp2, tol = tol, warmup = warmup, verbose = verbose)
  cp_list$cp2 <- cp2_dsn$parameters$cp
  mcp2 <- mean(exp(cp2_dsn$lp))

  ## one changepoint model
  cp1_dsn <- cp1_gibbs_v2(data = data, iter = iter, start.vals.left = start.vals$cp1$left, start.vals.right = start.vals$cp1$right,
                          prop_var_left = prop_var$cp1$left, prop_var_right = prop_var$cp1$right, cp_prop_var = cp_prop_var$cp1, tol = tol, warmup = warmup, verbose = verbose)
  cp_list$cp1 <- list("ppleft" = cp1_dsn$ppleft, "ppright" = cp1_dsn$ppright, "left" = cp1_dsn$left_parameters$cp, "right" = cp1_dsn$right_parameters$cp)
  mcp1 <- 0.5 * mean(exp(cp1_dsn$lp$left)) + 0.5 * mean(exp(cp1_dsn$lp$right))

  ## zero changepoint model
  cp0_dsn <- cp0_gibbs(data = data, iter = iter, start.vals = start.vals$cp0, prop_var = prop_var$cp0, tol = tol, warmup = warmup, verbose = verbose)
  mcp0 <- mean(cp0_dsn$lp)

  ## posterior cp probabilities
  ratio02 <- exp(log(prior_numcp[1]) + log(mcp0) - log(prior_numcp[3]) - log(mcp2))
  ratio12 <- exp(log(prior_numcp[2]) + log(mcp1) - log(prior_numcp[3]) - log(mcp2))

  p2 <- 1/(ratio02 + ratio12 + 1)
  p1 <- ratio12 * p2
  p0 <- ratio02 * p2

  post_numcp <- c(p0,p1,p2)

  avg_lp <- list("cp0" = mean(cp0_dsn$lp), "cp1" = 0.5 * mean(cp1_dsn$lp$left) + 0.5 * mean(cp1_dsn$lp$right), "cp2" = mean(cp2_dsn$lp))
  avg_lp_left <- mean(cp1_dsn$lp$left)
  avg_lp_right <- mean(cp1_dsn$lp$right)

  return(list("posterior_numcp" = post_numcp, "posterior_cp" = cp_list,
              "cp_mean" = list("2cp" = apply(X = cp_list$cp2, MARGIN = 2, FUN = mean),
                               "1cp" = list("left" = mean(cp_list$cp1$left), "right" = mean(cp_list$cp1$right))), "avg_lp" = avg_lp, "avg_lp_1cp" = c(avg_lp_left,avg_lp_right)))

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


