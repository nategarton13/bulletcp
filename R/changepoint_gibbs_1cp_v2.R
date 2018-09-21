cp1_gibbs_v2 <- function(data, iter, start.vals.left, start.vals.right, prop_var_left, prop_var_right, cp_prop_var, tol = 10, warmup = 5000, verbose = FALSE)
{
  ##data is a data frame with column x and column y

  ## run left cp algorithm first
  left_cp_out <- cp1_gibbs_left(data = data, iter = iter, warmup = warmup, start.vals = start.vals.left, prop_var = prop_var_left, cp_prop_var = cp_prop_var, verbose = verbose)

  ## run right cp algorithm
  right_cp_out <- cp1_gibbs_right(data = data, iter = iter, warmup = warmup, start.vals = start.vals.right, prop_var = prop_var_right, cp_prop_var = cp_prop_var, verbose = verbose)

  ## compute posterior probabilities of left or right changepoint
  # prior <- 0.5
  # pleft <- (prior * mean(left_cp_out$lp)) / (prior * mean(left_cp_out$lp) + (1-prior) * mean(right_cp_out$lp))
  # pright <- 1 - pleft
  pleft <- NA
  pright <- NA

  ## estimated changepoint
  est_left_cp <- mean(left_cp_out$parameters$cp)
  est_right_cp <- mean(right_cp_out$parameters$cp)

  return(list("left_parameters" = left_cp_out$parameters, "right_parameters" = right_cp_out$parameters,
              "accept" = list("left" = left_cp_out$accept, "right" = right_cp_out$accept),
              "lp" = list("left" = left_cp_out$lp, "right" = right_cp_out$lp),
              "ppleft" = pleft, "ppright" = pright, "estimated_cp" = list("left" = est_left_cp, "right" = est_right_cp))
         )
}

