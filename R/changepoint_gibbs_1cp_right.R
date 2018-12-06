#'  Estimate a posterior distribution of data conditional on a left groove and no right groove.
#'
#' This function runs a random walk metropolis within Gibbs algorithm to estimate the posterior distribution
#' of the value of the changepoint as well as the parameters fit in each multivariate normal distribution
#' on either side of the changepoint. The covariance matrices are both based on the exponential covariance function.
#' This functions assumes equally spaced locations ("x" values in the "data" argument). The distribution
#' to the right of the changepoint has a mean that is a linear function of the distance from the center of the data.
#' Note that this function is identical to the 1cp_left function, and more thorough
#' documentation is in that file.
#' @param data Data frame with columns "x" and "y." "x" is a column of the locations of the
#' observed residual values, y.
#' @param iter Number of interations after warmup.
#' @param start.vals List with elements "sigma", "l", "cp", "beta", and "intercept." "sigma" and "l"
#'   are 2 element vectors where the first element is for the data to the left of the changepoint.
#'   "cp" is the changepoint starting value. "beta" and "intercept" are the slope and intercept
#'   starting values for the mean of the data model to the left of the changepoint.
#'  which parameterize the covariance matrix.
#' @param prop_var A two element list of the proposal variance-covariance matrices for the random
#' walk metropolis algorithm(s). The first element is for the data to the left of the changepoint.
#' @param cp_prop_var The proposal variance for the changepoint.
#' @param tol_edge This parameter controls how close changepoint proposals can be to the edge of the data
#' before getting automatically rejected. For example, a value of 10 means that the changepoint will be
#' automatically rejected if the proposal is within a distance of 10 x-values from either edge.
#' @param warmup The number of initial iterations which serves two purposes: the first is to allow the
#' algorithm to wander to the area of most mass, and the second is to tune the proposal variance.
#' @param verbose Logical value indicating whether to print the iteration number and the parameter proposals.
#' @return A named list. "parameters" is a list of named parameter values each of which is a vector of length
#' "iter". "accept" gives the proportion of accepted proposals after warmup. "lp" is a vector of
#' values of the log data pdf at each sampled parameter value. "gp_prop_var" and "cp_prop_var" are
#' the tuned proposal variances for the metropolis steps.
#' @export
cp1_gibbs_right <- function(data, iter, start.vals, prop_var, cp_prop_var, tol_edge = 50, warmup = 500, verbose = FALSE)
{
  ##data is a data frame with column x and column y

  ## initialize parameter list
  par <- list()
  par$sigma <- matrix(nrow = warmup + 1, ncol = 2) ## the variance of the GP
  par$sigma[1,] <- start.vals$sigma

  par$l <- matrix(nrow = warmup + 1, ncol = 2) ## length scale of the GP
  par$l[1,] <- start.vals$l

  # par$tau <- matrix(nrow = warmup + 1, ncol = 2) ## nugget of the data model
  # par$tau[1,] <- start.vals$tau

  par$cp <- matrix(nrow = warmup + 1, ncol = 1) ## changepoint locations
  par$cp[1,] <- start.vals$cp

  par$beta <- matrix(nrow = warmup + 1, ncol = 1) ## regression coefficients for the GEAs
  par$beta[1,] <- start.vals$beta ## each row is the slope coefficient

  par$intercept <- matrix(nrow = warmup + 1, ncol = 1) ## regression intercept for the GEAs
  par$intercept[1,] <- start.vals$intercept ## each row is the slope coefficient

  ## range on the x-axis of data
  interval <- range(data$x)

  ## current values of parameters
  sigma <- start.vals$sigma
  l <- start.vals$l
  # tau <- start.vals$tau
  cp <- start.vals$cp
  beta <- start.vals$beta
  intercept <- start.vals$intercept


  ## initialize acceptance rates
  accept <- list()
  accept$gp_par <- matrix(data = c(0,0), nrow = 1, ncol = 2)
  accept$cp <- 0

  ## gibbs warmup iterations
  for(i in 1:(warmup))
  {
    xrange <- matrix(nrow = 2, ncol = 2)
    xrange[1,] <- c(interval[1], cp[1])
    xrange[2,] <- c(cp[1], interval[2])

    ## given changepoints make proposal for MH steps for GP parameters
    for(j in 1:2)
    {
      if(j == 2)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j], beta, intercept), sigma = prop_var[[j]]))
      }
      if(j == 1)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j]), sigma = prop_var[[j]]))
      }
      if(verbose == TRUE)
      {
        print(paste("iteration: ",i))
        print(paste(j,"-th GP parameter proposal: ", prop))
      }

      ## skip this chunk of data if the proposals result in values producing zero density
      if(j == 2)
      {
        if(any(prop[1:2] <= 0) || prop[3] <= 0)
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }
      if(j == 1)
      {
        if(any(prop <= 0))
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }

      temp_dat <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$y

      ## proposal doesn't appear because it should cancel
      if(j == 2)
      {
        med <- median(data$x)
        mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med)/(xrange[2,2] - xrange[1,1])) * beta[1] + intercept
        prop_mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med) / (xrange[2,2] - xrange[1,1])) * prop[3] + prop[4]

        log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = prop_mu, sigma = prop[1], l = prop[2]) + ## likelihood
          dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) + ## length scale
          dnorm(x = prop[3], mean = 0, sd = 10, log = TRUE) + ## slope
          dnorm(x = prop[4], mean = 0, sd = 10, log = TRUE) + ## intercept
          dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) - ## marginal standard deviation
          (lognormal_ou_pdf(x = temp_dat, mu = mu, sigma = sigma[j], l = l[j]) + ## likelihood
             dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) + ## length scale
             dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
             dnorm(x = beta, mean = 0, sd = 10, log = TRUE) + ## slope
             dnorm(x = intercept, mean = 0, sd = 10, log = TRUE)) ## intercept

        if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
        {
          sigma[j] <- prop[1]
          l[j] <- prop[2]
          # tau[j] <- prop[3]
          beta <- prop[3]
          intercept <- prop[4]
        }
      }
      if(j == 1)
      {
        log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2]) +
          dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) +
          dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) -
          (lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma[j], l = l[j]) +
             dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) +
             dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE))
        if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
        {
          sigma[j] <- prop[1]
          l[j] <- prop[2]
          # tau[j] <- prop[3]
        }
      }
    }
    ## update GP parameters
    par$sigma[i + 1,] <- sigma
    par$l[i + 1,] <- l
    # par$tau[i + 1,] <- tau
    par$beta[i + 1,] <- beta
    par$intercept[i + 1,] <- intercept

    ## sample from changepoint distribution given GP parameters
    prop <- as.numeric(rnorm(n = 1, mean = cp, sd = sqrt(cp_prop_var)))
    if(verbose == TRUE)
    {
      print(paste(i,"-th CP proposal: ", prop))
    }
    if(prop <= tol_edge + interval[1] || prop >= -tol_edge + interval[2])
    {
      par$cp[i + 1,] <- cp
    }
    else{
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y

      prop_temp_dat1 <- data[data$x <= prop & data$x > interval[1], ]$y
      prop_temp_dat2 <- data[data$x < interval[2] & data$x > prop, ]$y

      med2 <- median(data$x)
      mu2 <- ((data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$x - med2) / (xrange[2,2] - xrange[1,1])) * beta[1] + intercept
      mu1 <- rep(0, times = length(temp_dat1))

      prop_med2 <- median(data$x)
      prop_mu2 <- ((data[data$x > prop & data$x <= interval[2], ]$x - prop_med2) / (xrange[2,2] - xrange[1,1])) * beta[1] + intercept
      prop_mu1 <- rep(0, times = length(prop_temp_dat1))

      log_accept_ratio <- lognormal_ou_pdf(x = prop_temp_dat1, mu = prop_mu1, sigma = sigma[1], l = l[1]) +
        lognormal_ou_pdf(x = prop_temp_dat2, mu = prop_mu2, sigma = sigma[2], l = l[2]) -
        (lognormal_ou_pdf(x = temp_dat1, mu = mu1, sigma = sigma[1], l = l[1]) +
           lognormal_ou_pdf(x = temp_dat2, mu = mu2, sigma = sigma[2], l = l[2]))
      if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
      {
        cp <- prop
      }
    }
    par$cp[i + 1,] <- cp
    #print(i)
  }
  ###########################################################
  ## End warmup
  ###########################################################
  ## tune metropolis proposal variances
  prop_var[[2]] <- 2.4^2 * var(cbind(par$sigma[round(warmup/2):warmup,2], par$l[round(warmup/2):warmup,2], par$beta[round(warmup/2):warmup,1], par$intercept[round(warmup/2):warmup,1])) / 4 + 1e-1 * diag(4)
  prop_var[[1]] <- 2.4^2 * var(cbind(par$sigma[round(warmup/2):warmup,1], par$l[round(warmup/2):warmup,1])) / 2 + 1e-1 * diag(2)
  cp_prop_var <- 2.4^2 * var(par$cp[round(warmup/2):warmup,]) + 1

  ## reinitialize parameter list
  lp <- numeric()
  lpost <- numeric() ## log posterior values
  par <- list()
  par$sigma <- matrix(nrow = iter + 1, ncol = 2) ## the variance of the GP
  par$sigma[1,] <- sigma

  par$l <- matrix(nrow = iter + 1, ncol = 2) ## length scale of the GP
  par$l[1,] <- l

  # par$tau <- matrix(nrow = iter + 1, ncol = 2) ## nugget of the data model
  # par$tau[1,] <- tau

  par$beta <- matrix(nrow = iter + 1, ncol = 1) ## regression coefficient for the GEAs
  par$beta[1,] <- beta ## each row is the two slope coefficients

  par$intercept <- matrix(nrow = iter + 1, ncol = 1) ## regression intercept
  par$intercept[1,] <- intercept

  par$cp <- matrix(nrow = iter + 1, ncol = 1) ## changepoint locations
  par$cp[1,] <- cp

  ## gibbs sampling iterations
  for(i in 1:(iter))
  {
    xrange <- matrix(nrow = 2, ncol = 2)
    xrange[1,] <- c(interval[1], cp[1])
    xrange[2,] <- c(cp[1], interval[2])

    ## given changepoints make proposal for MH steps for GP parameters
    for(j in 1:2)
    {
      if(j == 2)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j], beta, intercept), sigma = prop_var[[j]]))
      }
      if(j == 1)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j]), sigma = prop_var[[j]]))
      }
      if(verbose == TRUE)
      {
        print(paste("iteration: ",i))
        print(paste(j,"-th GP parameter proposal: ", prop))
      }

      ## skip this chunk of data if the proposals result in values producing zero density
      if(j == 2)
      {
        if(any(prop[1:2] <= 0) || prop[3] <= 0)
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }
      if(j == 1)
      {
        if(any(prop <= 0))
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }

      temp_dat <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$y

      ## proposal doesn't appear because it should cancel
      if(j == 2)
      {
        med <- median(data$x)
        mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med)/(xrange[2,2] - xrange[1,1])) * beta[1] + intercept
        prop_mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med) / (xrange[2,2] - xrange[1,1])) * prop[3] + prop[4]

        log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = prop_mu, sigma = prop[1], l = prop[2]) + ## likelihood
          dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) + ## length scale
          dnorm(x = prop[3], mean = 0, sd = 10, log = TRUE) + ## slope
          dnorm(x = prop[4], mean = 0, sd = 10, log = TRUE) + ## intercept
          dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) - ## marginal standard deviation
          (lognormal_ou_pdf(x = temp_dat, mu = mu, sigma = sigma[j], l = l[j]) + ## likelihood
             dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) + ## length scale
             dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
             dnorm(x = beta, mean = 0, sd = 10, log = TRUE) + ## slope
             dnorm(x = intercept, mean = 0, sd = 10, log = TRUE)) ## intercept

        if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
        {
          accept$gp_par[1,j] <- accept$gp_par[1,j] + 1/iter
          sigma[j] <- prop[1]
          l[j] <- prop[2]
          # tau[j] <- prop[3]
          beta <- prop[3]
          intercept <- prop[4]
        }
      }
      if(j == 1)
      {
        log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = prop[1], l = prop[2]) +
          dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) +
          dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) -
          (lognormal_ou_pdf(x = temp_dat, mu = rep(0, times = length(temp_dat)), sigma = sigma[j], l = l[j]) +
             dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) +
             dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE))
        if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
        {
          accept$gp_par[1,j] <- accept$gp_par[1,j] + 1/iter
          sigma[j] <- prop[1]
          l[j] <- prop[2]
          # tau[j] <- prop[3]
        }
      }
    }
    ## update GP parameters
    par$sigma[i + 1,] <- sigma
    par$l[i + 1,] <- l
    # par$tau[i + 1,] <- tau
    par$beta[i + 1,] <- beta
    par$intercep[i + 1,] <- intercept

    ## sample from changepoint distribution given GP parameters
    prop <- as.numeric(rnorm(n = 1, mean = cp, sd = sqrt(cp_prop_var)))
    if(verbose == TRUE)
    {
      print(paste(i,"-th CP proposal: ", prop))
    }
    if(prop <= tol_edge + interval[1] || prop >= -tol_edge + interval[2])
    {
      par$cp[i + 1,] <- cp
      lp[i] <- (lognormal_ou_pdf(x = temp_dat1, mu = rep(0, times = length(temp_dat1)), sigma = sigma[1], l = l[1]) +
                  lognormal_ou_pdf(x = temp_dat2, mu = rep(0, times = length(temp_dat2)), sigma = sigma[2], l = l[2]))
      lpost[i] <- lp[i] + dgamma(x = l[1], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[1], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
        dnorm(x = beta, mean = 0, sd = 10, log = TRUE) + ## slope
        dnorm(x = intercept, mean = 0, sd = 10, log = TRUE) +
        dgamma(x = l[2], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[2], mean = 0, sd = 1, log = TRUE) ## marginal standard deviation
    }
    else{
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y

      prop_temp_dat1 <- data[data$x <= prop & data$x > interval[1], ]$y
      prop_temp_dat2 <- data[data$x < interval[2] & data$x > prop, ]$y

      med2 <- median(data$x)
      mu2 <- ((data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$x - med2) / (xrange[2,2] - xrange[1,1])) * beta[1] + intercept
      mu1 <- rep(0, times = length(temp_dat1))

      prop_med2 <- median(data$x)
      prop_mu2 <- ((data[data$x > prop & data$x <= interval[2], ]$x - prop_med2) / (xrange[2,2] - xrange[1,1])) * beta[1] + intercept
      prop_mu1 <- rep(0, times = length(prop_temp_dat1))

      log_accept_ratio <- lognormal_ou_pdf(x = prop_temp_dat1, mu = prop_mu1, sigma = sigma[1], l = l[1]) +
        lognormal_ou_pdf(x = prop_temp_dat2, mu = prop_mu2, sigma = sigma[2], l = l[2]) -
        (lognormal_ou_pdf(x = temp_dat1, mu = mu1, sigma = sigma[1], l = l[1]) +
           lognormal_ou_pdf(x = temp_dat2, mu = mu2, sigma = sigma[2], l = l[2]))

      lp[i] <- (lognormal_ou_pdf(x = temp_dat1, mu = rep(0, times = length(temp_dat1)), sigma = sigma[1], l = l[1]) +
                  lognormal_ou_pdf(x = temp_dat2, mu = rep(0, times = length(temp_dat2)), sigma = sigma[2], l = l[2]))

      lpost[i] <- lp[i] + dgamma(x = l[1], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[1], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
        dnorm(x = beta, mean = 0, sd = 10, log = TRUE) + ## slope
        dnorm(x = intercept, mean = 0, sd = 10, log = TRUE) +
        dgamma(x = l[2], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[2], mean = 0, sd = 1, log = TRUE)
      if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
      {
        cp <- prop
        accept$cp <- accept$cp + 1/iter
        lp[i] <- lognormal_ou_pdf(x = prop_temp_dat1, mu = rep(0, times = length(prop_temp_dat1)), sigma = sigma[1], l = l[1]) +
          lognormal_ou_pdf(x = prop_temp_dat2, mu = rep(0, times = length(prop_temp_dat2)), sigma = sigma[2], l = l[2])
        lpost[i] <- lp[i] + dgamma(x = l[1], shape = 3, rate = 5, log = TRUE) + ## length scale
          dnorm(x = sigma[1], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
          dnorm(x = beta, mean = 0, sd = 10, log = TRUE) + ## slope
          dnorm(x = intercept, mean = 0, sd = 10, log = TRUE) +
          dgamma(x = l[2], shape = 3, rate = 5, log = TRUE) + ## length scale
          dnorm(x = sigma[2], mean = 0, sd = 1, log = TRUE)
      }
    }
    par$cp[i + 1,] <- cp
    #print(i)
  }

  return(list("parameters" = par, "accept" = accept, "lp" = lp, "lpost" = lpost,"gp_prop_var" = prop_var, "cp_prop_var" = cp_prop_var))
}

# prop_var_right <- list(diag(c(0.2,0.2,0.2)), diag(c(0.2,0.2,0.2,0.2)))
# lopo <- cp1_gibbs_right(data = test.impute, start.vals = start.vals.right, prop_var = prop_var_right, cp_prop_var = 10^2, iter = 5000, warmup = 10000)
# lopo$accept
# plot(1:5001, lopo$parameters$cp)
# plot(1:5001, lopo$parameters$beta)
# plot(1:5001, lopo$parameters$l[,2])
# plot(1:5001, lopo$parameters$l[,1])
# plot(1:5001, lopo$parameters$sigma[,2])
