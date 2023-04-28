# logistic binomial model

# log likelihood
lbm_logl = function(beta1, beta2, gamma, Y, n, Z, X) {


  eta = X %*% gamma
  p = ilogit(eta)

  mu1 = ilogit(Z %*% (beta1 + beta2))
  mu2 = ilogit(Z %*% beta1)

  f1 = dbinom(Y, n, mu1)
  f2 = dbinom(Y, n, mu2)

  ll = sum(log(p * f1 + (1 - p) * f2))

  return(ll)
}

# generate a sample from the logistic binomial mixture model given mean parameters
# n: sample size
# mu1: Z * (beta1 + beta2)
# mu2: Z * beta1
# eta: X * gamma
# sizes: number of trials corresponding to mu1 and mu2
# Return: outcome vector Y
rlbm = function(n, mu1, mu2, eta, sizes) {

  # randomly generate binary class based on logistic regression
  class = rbinom(n, 1, ilogit(eta))

  # sample from mixture model with known class
  Y = numeric(n)
  for (i in 1:n) {
    if (class[i] == 1) {
      Y[i] = rbinom(1, sizes[i], ilogit(mu1[i]))
    }
    else if (class[i] == 0) {
      Y[i] = rbinom(1, sizes[i], ilogit(mu2[i]))
    }
  }

  # return outcome
  return(Y)

}

# fit the logistic binomiL model using metaheuristics
# Y: outcome vector
# Z: design matrix for mixture regression
# X: design matrix for logistic regression
# iter: maximum number of iterations for the algorithm
# swarm: population size for algorithm
# algorithm: one of the options from metaheuristicOpt
# Return: a named list of parameters
# This function assumes that the second column of Z contains the treatment variable
# therefore, we make the assumption that beta2[2] > 0 for identifibility
# note the sign difference from lnm
# also assume that there intercepts in both regression equations
fit_lbm = function(Y, sizes, Z, X, iter, swarm, algorithm) {

  # get dimensions
  n = length(Y)
  q1 = ncol(Z)
  q2 = ncol(X)

  # set up objective function
  obj_fun = function(param) {

    beta1 = param[1:q1]
    beta2 = param[(q1 + 1):(2*q1)]
    gamma = param[(2*q1+1):(2*q1 + q2)]
    LL = lbm_logl(beta1, beta2, gamma, Y, sizes, Z, X)

    # deal with missing
    if (is.na(LL))
      return(-Inf)
    else
      return(LL)
  }

  # set bounds for problem
  minY = min(Y)
  maxY = max(Y)
  rangeY = maxY - minY
  beta1_bounds = matrix(c(
    -10, 10, rep(c(-5, 5), q1 - 1)
  ), nrow = 2)
  beta2_bounds = matrix(c(
    -5, 5, 0, 5, rep(c(-5, 5), q1 - 2)
  ), nrow = 2)
  gamma_bounds = matrix(c(
    -10, 10, rep(c(-5, 5), q2 - 1)
  ), nrow = 2)

  bounds = cbind(beta1_bounds, beta2_bounds, gamma_bounds)

  # call metaheuristics library to maximize likelihood
  control = list()
  control$maxIter = iter
  control$numPopulation = swarm
  out = metaheuristicOpt::metaOpt(
    FUN = obj_fun,
    optimType = "MAX",
    algorithm = algorithm,
    numVar = 2*q1 + q2,
    rangeVar = bounds,
    control = control
  )


  # extract and return parameters
  beta1 = out$result[1:q1]
  beta2 = out$result[(q1 + 1):(2*q1)]
  gamma = out$result[(2*q1+1):(2*q1 + q2)]
  return(list(
    ll = out$optimumValue[1],
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma
  ))

}
