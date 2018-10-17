#'  Estimate a posterior distribution of data conditional that there are two grooves.
#'
#' This function runs a random walk metropolis within Gibbs algorithm to estimate the posterior distribution
#' of the value of the changepoints as well as the parameters fit in each multivariate normal distribution
#' on either side of each changepoint. The covariance matrices are based on the exponential covariance function.
#' This functions assumes equally spaced locations ("x" values in the "data" argument). The distribution
#' to the right of the right most changepoint and to the left of the left most changepoint have
#' means that are a linear function of the distance from the center of the data. The slope is
#' constrained to be negative in the left case and positive in the right case. The models fit
#' to the groove engraved areas are exactly the same as in the one changepoint case. Thus, this algorithm
#' only differs in that there are three segments of data to deal with as opposed to two.
#' @param data Data frame with columns "x" and "y." "x" is a column of the locations of the
#' observed residual values, y.
#' @param iter Number of interations after warmup.
#' @param start.vals Starting values for the changepoint algorithm.
#' List with elements "sigma", "l", "cp", "beta", and "intercept." "sigma" and "l"
#'   are 3 element vectors where the first element is for the data on the left groove.
#'   The second element is for the land engraved area, and the third element is for the right groove.
#'   "cp" is the vector of changepoint starting values. "beta" and "intercept" are two element vectors
#'   of the slope and intercept for the left and right groove engraved area respectively.
#' @param prop_var A three element list of the proposal variance-covariance matrices for the random
#' walk Metropolis algorithm(s). The first element is for the left groove engraved area.
#' The second element is for the land engraved area, and the third element is for the right engraved area.
#' @param cp_prop_var The proposal variance-covariance matrix for the changepoints.
#' @param tol_edge This parameter controls how close changepoint proposals can be to the edge of the data
#' before getting automatically rejected. For example, a value of 10 means that the changepoint will be
#' automatically rejected if either of the proposal changepoints is within a distance of 10 x-values from either edge.
#' @param tol_cp This parameter controls how close changepoint proposals can be to each other
#' before getting automatically rejected. For example, a value of 10 means that the changepoint will be
#' automatically rejected if either of the proposal changepoints is within a distance of 10 x-values from either each other.
#' @param warmup The number of initial iterations which serves two purposes: the first is to allow the
#' algorithm to wander to the area of most mass, and the second is to tune the proposal variance.
#' @param verbose Logical value indicating whether to print the iteration number and the parameter proposals.
#' @return A named list containing the sampled parameters, acceptance rates for the Metropolis steps,
#' log likelihood values, and proposal variance for the changepoints.
#' @importFrom mvtnorm rmvnorm
#' @export
cp2_gibbs_v2 <- function(data, iter, start.vals, prop_var, cp_prop_var, tol_edge = 50, tol_cp = 1000, warmup = 5000, verbose = FALSE)
{
  ##data is a data frame with column x and column y

  ## initialize parameter list
  par <- list()
  par$sigma <- matrix(nrow = warmup + 1, ncol = 3) ## the variance of the GP
  par$sigma[1,] <- start.vals$sigma

  par$l <- matrix(nrow = warmup + 1, ncol = 3) ## length scale of the GP
  par$l[1,] <- start.vals$l

  # par$tau <- matrix(nrow = warmup + 1, ncol = 3) ## nugget of the data model
  # par$tau[1,] <- start.vals$tau

  par$cp <- matrix(nrow = warmup + 1, ncol = 2) ## changepoint locations
  par$cp[1,] <- start.vals$cp

  par$beta <- matrix(nrow = warmup + 1, ncol = 2) ## regression coefficients for the GEAs
  par$beta[1,] <- start.vals$beta ## each row is the two slope coefficients

  par$intercept <- matrix(nrow = warmup + 1, ncol = 2) ## regression intercepts for the GEAs
  par$intercept[1,] <- start.vals$intercept ## each row is the two regression intercepts for GEAs

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
  accept$gp_par <- matrix(data = c(0,0,0), nrow = 1, ncol = 3)
  accept$cp <- 0

  ## gibbs warmup iterations
  for(i in 1:(warmup))
  {
    xrange <- matrix(nrow = 3, ncol = 2)
    xrange[1,] <- c(interval[1], cp[1])
    xrange[2,] <- c(cp[1], cp[2])
    xrange[3,] <- c(cp[2], interval[2])

    ## given changepoints make proposal for MH steps for GP parameters
    for(j in 1:3)
    {
      if(j == 1)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j], beta[1], intercept[1]), sigma = prop_var[[j]]))
      }
      if(j == 3)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j], beta[2], intercept[2]), sigma = prop_var[[j]]))
      }
      if(j == 2)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j]), sigma = prop_var[[j]]))
      }
      if(verbose == TRUE)
      {
        print(paste("iteration: ",i))
        print(paste(j,"-th GP parameter proposal: ", prop))
      }

      ## skip this chunk of data if the proposals result in values producing zero density
      if(j == 1)
      {
        if(any(prop[1:2] <= 0) || prop[3] >= 0)
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }
      if(j == 2)
      {
        if(any(prop <= 0))
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }
      if(j == 3)
      {
        if(any(prop[1:2] <= 0) || prop[3] <= 0)
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }

        temp_dat <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$y

        ## proposal doesn't appear because it should cancel
        if(j == 1)
        {
          med <- median(data$x)
          mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med)/(xrange[3,2] - xrange[1,1])) * beta[1] + intercept[1]
          prop_mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x) / (xrange[3,2] - xrange[1,1])) * prop[3] + prop[4]

          log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = prop_mu, sigma = prop[1], l = prop[2]) + ## likelihood
            dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) + ## length scale
            dnorm(x = prop[3], mean = 0, sd = 10, log = TRUE) + ## slope
            dnorm(x = prop[4], mean = 0, sd = 10, log = TRUE) + ## intercept
            dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) - ## marginal standard deviation
            (lognormal_ou_pdf(x = temp_dat, mu = mu, sigma = sigma[j], l = l[j]) + ## likelihood
               dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) + ## length scale
               dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
               dnorm(x = beta[1], mean = 0, sd = 10, log = TRUE) + ## intercept
               dnorm(x = intercept[1], mean = 0, sd = 10, log = TRUE)) ## slope

          if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
          {
            sigma[j] <- prop[1]
            l[j] <- prop[2]
            # tau[j] <- prop[3]
            beta[1] <- prop[3]
            intercept[1] <- prop[4]
          }
        }
        if(j == 3)
        {
          med <- median(data$x)
          mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med)/(xrange[3,2] - xrange[1,1])) * beta[2] + intercept[2]
          prop_mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x) / (xrange[3,2] - xrange[1,1])) * prop[3] + prop[4]

          log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = prop_mu, sigma = prop[1], l = prop[2]) + ## likelihood
            dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) + ## length scale
            dnorm(x = prop[3], mean = 0, sd = 10, log = TRUE) + ## slope
            dnorm(x = prop[4], mean = 0, sd = 10, log = TRUE) + ## intercept
            dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) - ## marginal standard devivation
            (lognormal_ou_pdf(x = temp_dat, mu = mu, sigma = sigma[j], l = l[j]) + ## likelihood
               dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) + ## length scale
               dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
               dnorm(x = intercept[2], mean = 0, sd = 10, log = TRUE) + ## intercept
               dnorm(x = beta[2], mean = 0, sd = 10, log = TRUE)) ## slope

          if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
          {
            sigma[j] <- prop[1]
            l[j] <- prop[2]
            # tau[j] <- prop[3]
            beta[2] <- prop[3]
            intercept[2] <- prop[4]
          }
        }
        if(j == 2)
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
    prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = cp, sigma = cp_prop_var))
    if(verbose == TRUE)
    {
      print(paste(i,"-th CP proposal: ", prop))
    }
    if(prop[1] >= prop[2] - tol_cp || prop[1] <= tol_edge + interval[1] || prop[2] >= -tol_edge + interval[2])
    {
      par$cp[i + 1,] <- cp
    }
    else{
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y
      temp_dat3 <- data[data$x <= xrange[3,2] & data$x > xrange[3,1], ]$y

      prop_temp_dat1 <- data[data$x <= prop[1] & data$x > interval[1], ]$y
      prop_temp_dat2 <- data[data$x <= prop[2] & data$x > prop[1], ]$y
      prop_temp_dat3 <- data[data$x <= interval[2] & data$x > prop[2], ]$y

      med1 <- median(data$x)
      med3 <- median(data$x)
      mu3 <- ((data[data$x <= xrange[3,2] & data$x > xrange[3,1], ]$x - med3) / (xrange[3,2] - xrange[1,1])) * beta[2] + intercept[2]
      mu1 <- ((data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$x - med1) / (xrange[3,2] - xrange[1,1])) * beta[1] + intercept[1]
      mu2 <- rep(0, times = length(temp_dat2))

      prop_med1 <- median(data$x)
      prop_med3 <- median(data$x)
      prop_mu3 <- ((data[data$x <= interval[2] & data$x > prop[2], ]$x - prop_med3) / (xrange[3,2] - xrange[1,1])) * beta[2] + intercept[2]
      prop_mu1 <- ((data[data$x <= prop[1] & data$x > interval[1], ]$x - prop_med1) / (xrange[3,2] - xrange[1,1])) * beta[1] + intercept[1]
      prop_mu2 <- rep(0, times = length(prop_temp_dat2))

      log_accept_ratio <- lognormal_ou_pdf(x = prop_temp_dat1, mu = prop_mu1, sigma = sigma[1], l = l[1]) +
        lognormal_ou_pdf(x = prop_temp_dat2, mu = prop_mu2, sigma = sigma[2], l = l[2]) +
        lognormal_ou_pdf(x = prop_temp_dat3, mu = prop_mu3, sigma = sigma[3], l = l[3]) -
        (lognormal_ou_pdf(x = temp_dat1, mu = mu1, sigma = sigma[1], l = l[1]) +
           lognormal_ou_pdf(x = temp_dat2, mu = mu2, sigma = sigma[2], l = l[2]) +
           lognormal_ou_pdf(x = temp_dat3, mu = mu3, sigma = sigma[3], l = l[3]))
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
  prop_var[[1]] <- 2.4^2 * var(cbind(par$sigma[round(warmup/2):warmup,1], par$l[round(warmup/2):warmup,1], par$beta[round(warmup/2):warmup,1], par$intercept[round(warmup/2):warmup,1])) / 4 + 1e-5 * diag(4)
  prop_var[[2]] <- 2.4^2 * var(cbind(par$sigma[round(warmup/2),2], par$l[round(warmup/2):warmup,2])) / 2 + 1e-5 * diag(2)
  prop_var[[3]] <- 2.4^2 * var(cbind(par$sigma[round(warmup/2):warmup,3], par$l[round(warmup/2):warmup,3], par$beta[round(warmup/2):warmup,2], par$intercept[round(warmup/2) : warmup, 2])) / 4 + 1e-5 * diag(4)
  cp_prop_var <- 2.4^2 * var(par$cp[round(warmup/2):warmup,]) / 2 + 1e-5 * diag(2)

  ## reinitialize parameter list
  lp <- numeric() ## the log likelihood
  lpost <- numeric() ## the log posterior values
  par <- list()
  par$sigma <- matrix(nrow = iter + 1, ncol = 3) ## the variance of the GP
  par$sigma[1,] <- sigma

  par$l <- matrix(nrow = iter + 1, ncol = 3) ## length scale of the GP
  par$l[1,] <- l

  # par$tau <- matrix(nrow = iter + 1, ncol = 3) ## nugget of the data model
  # par$tau[1,] <- tau

  par$cp <- matrix(nrow = iter + 1, ncol = 2) ## changepoint locations
  par$cp[1,] <- cp

  par$beta <- matrix(nrow = iter + 1, ncol = 2) ## slopes
  par$beta[1,] <- beta

  par$intercept <- matrix(nrow = iter + 1, ncol = 2) ## intercepts
  par$intercept[1,] <- intercept

  ## gibbs sampling iterations
  for(i in 1:(iter))
  {
    xrange <- matrix(nrow = 3, ncol = 2)
    xrange[1,] <- c(interval[1], cp[1])
    xrange[2,] <- c(cp[1], cp[2])
    xrange[3,] <- c(cp[2], interval[2])

    ## given changepoints make proposal for MH steps for GP parameters
    for(j in 1:3)
    {
      if(j == 1)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j], beta[1], intercept[1]), sigma = prop_var[[j]]))
      }
      if(j == 3)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j], beta[2], intercept[2]), sigma = prop_var[[j]]))
      }
      if(j == 2)
      {
        prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = c(sigma[j], l[j]), sigma = prop_var[[j]]))
      }
      if(verbose == TRUE)
      {
        print(paste("iteration: ",i))
        print(paste(j,"-th GP parameter proposal: ", prop))
      }

      ## skip this chunk of data if the proposals result in values producing zero density
      if(j == 1)
      {
        if(any(prop[1:2] <= 0) || prop[3] >= 0)
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }
      if(j == 2)
      {
        if(any(prop <= 0))
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }
      if(j == 3)
      {
        if(any(prop[1:2] <= 0) || prop[3] <= 0)
        {
          # par$sigma[i,j] <- sigma[j]
          # par$l[i,j] <- l[j]
          # par$tau[i,j] <- tau[j]
          next
        }
      }

        temp_dat <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$y

        ## proposal doesn't appear because it should cancel
        if(j == 1)
        {
          med <- median(data$x)
          mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med) / (xrange[3,2] - xrange[1,1])) * beta[1] + intercept[1]
          prop_mu <- ((data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x) / (xrange[3,2] - xrange[1,1])) * prop[3] + prop[4]

          log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = prop_mu, sigma = prop[1], l = prop[2]) + ## likelihood
            dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) + ## length scale
            dnorm(x = prop[3], mean = 0, sd = 10, log = TRUE) + ## slope
            dnorm(x = prop[4], mean = 0, sd = 10, log = TRUE) + ## intercept
            dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) - ## marginal standard deviation
            (lognormal_ou_pdf(x = temp_dat, mu = mu, sigma = sigma[j], l = l[j]) + ## likelihood
               dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) + ## length scale
               dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
               dnorm(x = intercept[1], mean = 0, sd = 10, log = TRUE) + ## intercept
               dnorm(x = beta[1], mean = 0, sd = 10, log = TRUE)) ## slope

          if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
          {
            accept$gp_par[1,j] <- accept$gp_par[1,j] + 1/iter
            sigma[j] <- prop[1]
            l[j] <- prop[2]
            # tau[j] <- prop[3]
            beta[1] <- prop[3]
            intercept[1] <- prop[4]
          }
        }
        if(j == 3)
        {
          med <- median(data$x)
          mu <- (data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x - med) / (xrange[3,2] - xrange[1,1]) * beta[2] + intercept[2]
          prop_mu <- data[data$x <= xrange[j,2] & data$x > xrange[j,1], ]$x / (xrange[3,2] - xrange[1,1]) * prop[3] + prop[4]

          log_accept_ratio <- lognormal_ou_pdf(x = temp_dat, mu = prop_mu, sigma = prop[1], l = prop[2]) + ## likelihood
            dgamma(x = prop[2], shape = 3, rate = 5, log = TRUE) + ## length scale
            dnorm(x = prop[3], mean = 0, sd = 10, log = TRUE) + ## slope
            dnorm(x = prop[4], mean = 0, sd = 10, log = TRUE) + ## intercept
            dnorm(x = prop[1], mean = 0, sd = 1, log = TRUE) - ## marginal standard devivation
            (lognormal_ou_pdf(x = temp_dat, mu = mu, sigma = sigma[j], l = l[j]) + ## likelihood
               dgamma(x = l[j], shape = 3, rate = 5, log = TRUE) + ## length scale
               dnorm(x = sigma[j], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
               dnorm(x = intercept[2], mean = 0, sd = 10, log = TRUE) + ## intercept
               dnorm(x = beta[2], mean = 0, sd = 10, log = TRUE)) ## slope

          if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
          {
            accept$gp_par[1,j] <- accept$gp_par[1,j] + 1/iter
            sigma[j] <- prop[1]
            l[j] <- prop[2]
            # tau[j] <- prop[3]
            beta[2] <- prop[3]
            intercept[2] <- prop[4]
          }
        }
        if(j == 2)
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
    par$intercept[i + 1,] <- intercept

    ## sample from changepoint distribution given GP parameters
    prop <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = cp, sigma = cp_prop_var))
    if(verbose == TRUE)
    {
      print(paste(i,"-th CP proposal: ", prop))
    }
    if(prop[1] >= prop[2] - tol_cp || prop[1] <= tol_edge + interval[1] || prop[2] >= -tol_edge + interval[2])
    {
      par$cp[i + 1,] <- cp
      lp[i] <- (lognormal_ou_pdf(x = temp_dat1, mu = mu1, sigma = sigma[1], l = l[1]) +
                  lognormal_ou_pdf(x = temp_dat2, mu = mu2, sigma = sigma[2], l = l[2]) +
                  lognormal_ou_pdf(x = temp_dat3, mu = mu3, sigma = sigma[3], l = l[3]))
      lpost[i] <- lp[i] + dgamma(x = l[1], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[1], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
        dnorm(x = intercept[1], mean = 0, sd = 10, log = TRUE) + ## intercept
        dnorm(x = beta[1], mean = 0, sd = 10, log = TRUE) +
        dgamma(x = l[2], shape = 3, rate = 5, log = TRUE) +
        dnorm(x = sigma[2], mean = 0, sd = 1, log = TRUE) +
        dgamma(x = l[3], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[3], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
        dnorm(x = intercept[2], mean = 0, sd = 10, log = TRUE) + ## intercept
        dnorm(x = beta[2], mean = 0, sd = 10, log = TRUE)

    }
    else{
      temp_dat1 <- data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$y
      temp_dat2 <- data[data$x <= xrange[2,2] & data$x > xrange[2,1], ]$y
      temp_dat3 <- data[data$x <= xrange[3,2] & data$x > xrange[3,1], ]$y

      prop_temp_dat1 <- data[data$x <= prop[1] & data$x > interval[1], ]$y
      prop_temp_dat2 <- data[data$x <= prop[2] & data$x > prop[1], ]$y
      prop_temp_dat3 <- data[data$x <= interval[2] & data$x > prop[2], ]$y

      med3 <- median(data$x)
      med1 <- median(data$x)
      mu3 <- (data[data$x <= xrange[3,2] & data$x > xrange[3,1], ]$x - med3) / (xrange[3,2] - xrange[1,1]) * beta[2] + intercept[2]
      mu1 <- (data[data$x <= xrange[1,2] & data$x > xrange[1,1], ]$x - med1) / (xrange[3,2] - xrange[1,1]) * beta[1] + intercept[1]
      mu2 <- rep(0, times = length(temp_dat2))

      prop_med3 <- median(data$x)
      prop_med1 <- median(data$x)
      prop_mu3 <- (data[data$x <= interval[2] & data$x > prop[2], ]$x - prop_med3) / (xrange[3,2] - xrange[1,1]) * beta[2] + intercept[2]
      prop_mu1 <- (data[data$x <= prop[1] & data$x > interval[1], ]$x - prop_med1) / (xrange[3,2] - xrange[1,1]) * beta[1] + intercept[1]
      prop_mu2 <- rep(0, times = length(prop_temp_dat2))

      log_accept_ratio <- lognormal_ou_pdf(x = prop_temp_dat1, mu = prop_mu1, sigma = sigma[1], l = l[1]) +
        lognormal_ou_pdf(x = prop_temp_dat2, mu = prop_mu2, sigma = sigma[2], l = l[2]) +
        lognormal_ou_pdf(x = prop_temp_dat3, mu = prop_mu3, sigma = sigma[3], l = l[3]) -
        (lognormal_ou_pdf(x = temp_dat1, mu = mu1, sigma = sigma[1], l = l[1]) +
           lognormal_ou_pdf(x = temp_dat2, mu = mu2, sigma = sigma[2], l = l[2]) +
           lognormal_ou_pdf(x = temp_dat3, mu = mu3, sigma = sigma[3], l = l[3]))

      lp[i] <- (lognormal_ou_pdf(x = temp_dat1, mu = mu1, sigma = sigma[1], l = l[1]) +
          lognormal_ou_pdf(x = temp_dat2, mu = mu2, sigma = sigma[2], l = l[2]) +
          lognormal_ou_pdf(x = temp_dat3, mu = mu3, sigma = sigma[3], l = l[3]))
      lpost[i] <- lp[i] + dgamma(x = l[1], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[1], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
        dnorm(x = intercept[1], mean = 0, sd = 10, log = TRUE) + ## intercept
        dnorm(x = beta[1], mean = 0, sd = 10, log = TRUE) +
        dgamma(x = l[2], shape = 3, rate = 5, log = TRUE) +
        dnorm(x = sigma[2], mean = 0, sd = 1, log = TRUE) +
        dgamma(x = l[3], shape = 3, rate = 5, log = TRUE) + ## length scale
        dnorm(x = sigma[3], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
        dnorm(x = intercept[2], mean = 0, sd = 10, log = TRUE) + ## intercept
        dnorm(x = beta[2], mean = 0, sd = 10, log = TRUE)

      if(log(runif(n = 1, min = 0, max = 1)) <= log_accept_ratio)
      {
        cp <- prop
        accept$cp <- accept$cp + 1/iter
        lp[i] <- lognormal_ou_pdf(x = prop_temp_dat1, mu = prop_mu1, sigma = sigma[1], l = l[1]) +
          lognormal_ou_pdf(x = prop_temp_dat2, mu = prop_mu2, sigma = sigma[2], l = l[2]) +
          lognormal_ou_pdf(x = prop_temp_dat3, mu = prop_mu3, sigma = sigma[3], l = l[3])
        lpost[i] <- lp[i] + dgamma(x = l[1], shape = 3, rate = 5, log = TRUE) + ## length scale
          dnorm(x = sigma[1], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
          dnorm(x = intercept[1], mean = 0, sd = 10, log = TRUE) + ## intercept
          dnorm(x = beta[1], mean = 0, sd = 10, log = TRUE) +
          dgamma(x = l[2], shape = 3, rate = 5, log = TRUE) +
          dnorm(x = sigma[2], mean = 0, sd = 1, log = TRUE) +
          dgamma(x = l[3], shape = 3, rate = 5, log = TRUE) + ## length scale
          dnorm(x = sigma[3], mean = 0, sd = 1, log = TRUE) + ## marginal standard deviation
          dnorm(x = intercept[2], mean = 0, sd = 10, log = TRUE) + ## intercept
          dnorm(x = beta[2], mean = 0, sd = 10, log = TRUE)
      }
    }
    par$cp[i + 1,] <- cp
    # print(lp[i])
    #print(i)
  }

  return(list("parameters" = par, "accept" = accept, "cp_prop_var" = cp_prop_var, "lp" = lp, "lpost" = lpost))
}
