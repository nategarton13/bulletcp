#'  Impute missing data.
#'
#' This function performs maximum likelihood estimation to estimate the variance parameters
#' in a Gaussian process with a squared exponential covariance function. These parameters could then be used
#' in the Gaussian process used for imputation.
#' @param y Numeric y vector of response values.
#' @param x Numeric x vector of locations used for the covariance function.
#' @param tol Tolerance level for the maximum likelihood procedure to fit the Gaussian process.
#' @return Standard optim output. The first optimized parameter value is the standard deviation
#'   the second is the length scale.
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats optim
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
#' fake_groove <- fake_groove[sample.int(n = nrow(fake_groove),
#'     size = round(0.8 * nrow(fake_groove)),
#'     replace = FALSE),]
#' plot(fake_groove$x, fake_groove$y)
#'
#'
#' # estimate the MLE's
#' mles <- mlgp(y = fake_groove$y, x = fake_groove$x)
#' @export
mlgp <- function(y,x,tol = 1e-6)
{
  ## squared exponential covariance function
  cov.fun <- function(x1, x2, par)
  {
    sigma <- exp(par[1])
    tau <- exp(par[2])
    return(sigma^2 * exp(-1/(2*tau^2) * (x1 - x2)^2))
  }

  ## function to make a covariance matrix based on a given covariance function cov.fun,
  ## a "nugget" to be used only if some values to be predicted occur at the same values as
  ## the data. x are the locations of the observed data, xpred are the locations of the
  ## data to be predicted. par is list("sigma", "l").
  make.cov.mat <- function(x, xpred = numeric(), cov.fun, par, eps = 0)
  {
    xfull <- c(x,xpred)
    temp <- matrix(nrow = length(xfull), ncol = length(xfull))
    for(i in 1:length(xfull))
    {
      for(j in 1:length(xfull))
      {
        temp[i,j] <- cov.fun(xfull[i], xfull[j], par = par) + 1*(i==j)*eps
      }
    }
    return(temp)
  }

  ## this is the multivariate normal objective function to do maximum likelihood with
  obj_fun <- function(par, ...)
  {
    ## make covariance matrix
    K <- make.cov.mat(x = x, xpred = numeric(), cov.fun = cov.fun, par = par, eps = 0)
    return(-dmvnorm(x = y, sigma = K, log = TRUE))
  }

  mle <- optim(par = c(0,1), fn = obj_fun, x = x, y = y)

  ## returns the values returned by optim
  ## par[1] is sigma, par[2] is the length scale
  return(mle)
}

#'  Impute missing data.
#'
#' This function imputes missing data based on a Gaussian process regression
#' @param y Numeric y vector of response values.
#' @param x Numeric x vector of locations used for the covariance function.
#' @param sigma Marginal standard deviation in the Gaussian process.
#' @param l Length scale parameter in the Gaussian process.
#' @return A data frame with columns "x" and "y" which contain the combined observed and imputed data.
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
#' fake_groove <- fake_groove[sample.int(n = nrow(fake_groove),
#'     size = round(0.8 * nrow(fake_groove)),
#'     replace = FALSE),]
#' fake_groove <- fake_groove[order(fake_groove$x),]
#' plot(fake_groove$x, fake_groove$y)
#'
#' # add NA values where the data are missing
#' x_na <- seq(from = min(fake_groove$x), to = max(fake_groove$x), by = min(fake_groove$x[2:nrow(fake_groove)] - fake_groove$x[1:(nrow(fake_groove) - 1)]))
#' x_na <- x_na[!round(x_na, digits = 2) %in% round(fake_groove$x, digits = 2)]
#' y_na <- rep(NA, times = length(x_na))
#' d_na <- data.frame("x" = x_na, "y" = y_na)
#' fake_groove <- rbind(fake_groove, d_na)
#' fake_groove <- fake_groove[order(fake_groove$x),]
#'
#' ## impute the data
#' full_data <- imputeGP(y = fake_groove$y, x = fake_groove$x, sigma = 0.9, l = 15)
#' head(full_data)
#' plot(full_data$x, full_data$y)
#' @export

imputeGP <- function(y,x, sigma, l)
{
  ## return a warning if there are no NA values in y
  if(any(is.na(y)) == FALSE)
  {
    print("Warning: no NA values detected in y. Returning initial data.")
    return(data.frame("x" = x, "y" = y))
  }

  ## squared exponential covariance function
  cov.fun <- function(x1, x2, par)
  {
    sigma <- par$sigma
    tau <- par$tau
    return(sigma^2 * exp(-1/(2*tau^2) * (x1 - x2)^2))
  }

  ## function to create covariance matrix
  make.cov.mat <- function(x, xpred, cov.fun, par, eps = 0)
  {
    xfull <- c(x,xpred)
    temp <- matrix(nrow = length(xfull), ncol = length(xfull))
    for(i in 1:length(xfull))
    {
      for(j in 1:length(xfull))
      {
        temp[i,j] <- cov.fun(xfull[i], xfull[j], par = par) + 1*(i==j)*eps
      }
    }
    return(temp)
  }

  ## function to evaluate the conditional mean of unobserved values at xpred
  ## based on data, y, at locations, x. mu is the mean of the entire vector and sigma
  ## is the covariance matrix of the entire vector.
  normal.cond.mean <- function(y, x, xpred, mu = rep(0, times = length(c(x,xpred))), sigma)
  {
    sigma21 <- sigma[(length(x) + 1):length(mu),1:length(x)]
    sigma11 <- sigma[1:length(x), 1:length(x)]
    return(mu[(length(x) + 1):length(mu)] + sigma21 %*% qr.solve(sigma11) %*% (y - mu[1:length(x)]))
  }

  ## function to evaluate the conditional variance of unobserved values at xpred
  ## based on data, y, at locations, x. mu is the mean of the entire vector and sigma
  ## is the covariance matrix of the entire vector.
  normal.cond.var <- function(y,x, xpred, sigma, mu)
  {
    sigma21 <- sigma[(length(x) + 1):length(mu),1:length(x)]
    sigma11 <- sigma[1:length(x), 1:length(x)]
    sigma22 <- sigma[(length(x) + 1):length(mu),(length(x) + 1):length(mu)]
    return(sigma22 - sigma21 %*% solve(sigma11) %*% t(sigma21))
  }

  ## split data into NA and non-NA data
  #y.NA <- y[is.na(y)]
  x.NA <- x[is.na(y)]
  y.ok <- y[!is.na(y)]
  x.ok <- x[!is.na(y)]
  x.ok <- x.ok[seq(from = 1, to = length(x.ok), by = 5)]
  y.ok <- y.ok[seq(from = 1, to = length(y.ok), by = 5)]

  ## make covariance matrix
  K <- make.cov.mat(x = x.ok, xpred = x.NA, cov.fun = cov.fun, par = list("sigma" = sigma, "tau" = l), eps = 1e-3)

  ## make predictions and return the data frame of observed and imputed values
  pred.y <- normal.cond.mean(y = y.ok, x = x.ok, xpred = x.NA, sigma = K)
  pred.df <- data.frame("x" = x, "y" = y)
  pred.df$y[is.na(y)] <- pred.y
  return(pred.df)
}

