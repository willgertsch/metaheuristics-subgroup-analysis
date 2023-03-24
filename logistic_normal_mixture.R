# logistic normal mixture model

ilogit = function(eta) {
  1/(1 + exp(-eta))
}

# log-likelihood
lnm_logl = function(beta1, beta2, sigma, gamma, Y, Z, X) {


  eta = X %*% gamma
  p = ilogit(eta)

  mu1 = Z %*% (beta1 + beta2)
  mu2 = Z %*% beta1

  f1 = dnorm(Y, mu1, sigma)
  f2 = dnorm(Y, mu2, sigma)

  ll = sum(log(p * f1 + (1 - p) * f2))

  return(ll)
}

# generate a sample from the logistic normal mixture model given mean parameters
# n: sample size
# mu1: Z * (beta1 + beta2)
# mu2: Z * beta1
# eta: X * gamma
# Return: outcome vector Y
rlnm = function(n, mu1, mu2, eta, sigma) {

  # randomly generate binary class based on logistic regression
  class = rbinom(n, 1, ilogit(eta))

  # sample from mixture model with known class
  Y = numeric(n)
  for (i in 1:n) {
    if (class[i] == 1) {
      Y[i] = rnorm(1, mu1[i], sigma)
    }
    else if (class[i] == 0) {
      Y[i] = rnorm(1, mu2[i], sigma)
    }
  }

  # return outcome
  return(Y)

}

# update the value of gamma during M step using NR
# X: matrix of covariates
# W: vector of posterior probabilities of class membership
# gamma_init: initial gamma values
# iter: number of iterations for RN
# Returns: updated gamma
max_Q1 = function(X, W, gamma_init, iter) {

  gamma = gamma_init # initial value
  for (i in 1:iter) {
    P = ilogit(X %*% gamma)
    grad = t(X) %*% (W - P)
    Omega = diag(c(P * (1 - P)))
    IM = t(X) %*% Omega %*% X
    gamma = gamma + solve(IM) %*% grad
    cat("Iteration: ", i,  ",  Gamma=", gamma, "\n")
  }
  return(gamma)
}

# update beta parameters and sigma using WLS steps
# Y: outcome vector
# Z: matrix of covariates
# W: Vector of posterior probabilities of class membership
# *_init: initial values for parameters
# Return: updated beta1, beta2, sigma
max_Q2 = function(Y, Z, W, beta1_init, beta2_init, sigma_init) {

  tbeta1_init = beta1_init + beta2_init # reparameterize
  tbeta2_init = beta1_init

  # beta estimation using WLS-like step
  # update tbeta1 and then tbeta2
  f1 = exp(-1/(2*sigma_init^2) * (Y - Z %*% tbeta1_init)^2)
  f2 = exp(-1/(2*sigma_init^2) * (Y - Z %*% tbeta2_init)^2)
  omega1 = 1/(W * f1 + (1 - W) * f2) * f1
  Omega1 = diag(c(omega1))
  hat_tbeta1 = solve(t(Z) %*% Omega1 %*% Z) %*% t(Z) %*% Omega1 %*% Y

  f1 = exp(-1/(2*sigma_init^2) * (Y - Z %*% hat_tbeta1)^2)
  omega2 = 1/(W * f1 + (1 - W) * f2) * f2
  Omega2 = diag(c(omega2))
  hat_tbeta2 = solve(t(Z) %*% Omega2 %*% Z) %*% t(Z) %*% Omega2 %*% Y

  # sigma estimation using WLS-like step
  # shortcut: averaging individual variance estimates
  sigma2_1 = sum(omega1*(Y - Z%*%hat_tbeta1)^2)/sum(omega1)
  sigma2_2 = sum(omega2*(Y - Z%*%hat_tbeta2)^2)/sum(omega2)
  hat_sigma = sqrt((sigma2_1 + sigma2_2)/2)

  return(list(beta1 = hat_tbeta2, beta2 = hat_tbeta1 - hat_tbeta2, sigma = hat_sigma))

}

# perform one step of the EM update
# returns list of updated parameter values
update_EM = function(Y, Z, X, beta1, beta2, sigma, gamma) {

  eta = X %*% gamma
  p = ilogit(eta)

  # E step
  mu1 = Z %*% (beta1 + beta2)
  mu2 = Z %*% beta1
  f1 = dnorm(Y, mu1, sigma)
  f2 = dnorm(Y, mu2, sigma)
  denoms = f1 + f2
  W = p * f1 / denoms

  # M step ####
  # gamma update
  # how many iterations of IRLS is good?
  gamma_t = max_Q1(X, W, gamma, 10)
  # update beta and sigma
  result = max_Q2(Y, Z, W, beta1, beta2, sigma)

  # return update values
  return(list(
    beta1 = result$beta1,
    beta2 = result$beta2,
    gamma = gamma_t,
    sigma = result$sigma
  ))

}

# fit the logistic normal model using metaheuristics
# Y: outcome vector
# Z: design matrix for mixture regression
# X: design matrix for logistic regression
# iter: maximum number of iterations for the algorithm
# swarm: population size for algorithm
# algorithm: one of the options from metaheuristicOpt
# Return: a named list of parameters
# This function assumes that the second column of Z contains the treatment variable
# therefore, we make the assumption that beta2[2] > 0 for identifibility
# also assume that there intercepts in both regression equations
fit_lnm = function(Y, Z, X, iter, swarm, algorithm) {

  # get dimensions
  n = length(Y)
  q1 = ncol(Z)
  q2 = ncol(X)

  # set up objective function
  obj_fun = function(param) {

    beta1 = param[1:q1]
    beta2 = param[(q1 + 1):(2*q1)]
    gamma = param[(2*q1+1):(2*q1 + q2)]
    sigma = param[2*q1 + q2 + 1]
    LL = lnm_logl(beta1, beta2, sigma, gamma, Y, Z, X)

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
    minY, maxY, rep(c(-rangeY, rangeY), q1 - 1)
  ), nrow = 2)
  beta2_bounds = matrix(c(
    -rangeY, rangeY, 0, rangeY, rep(c(-rangeY, rangeY), q1 - 2)
  ), nrow = 2)
  gamma_bounds = matrix(c(
    -10, 10, rep(c(-5, 5), q2 - 1)
  ), nrow = 2)
  sigma_bounds = matrix(c(
    0.1, sd(Y)
  ), nrow = 2)

  bounds = cbind(beta1_bounds, beta2_bounds, gamma_bounds, sigma_bounds)

  # call metaheuristics library to maximize likelihood
  control = list()
  control$maxIter = iter
  control$numPopulation = swarm
  out = metaheuristicOpt::metaOpt(
    FUN = obj_fun,
    optimType = "MAX",
    algorithm = algorithm,
    numVar = 2*q1 + q2 + 1,
    rangeVar = bounds,
    control = control
  )


  # extract and return parameters
  beta1 = out$result[1:q1]
  beta2 = out$result[(q1 + 1):(2*q1)]
  gamma = out$result[(2*q1+1):(2*q1 + q2)]
  sigma = out$result[2*q1 + q2 + 1]
  return(list(
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma,
    sigma = sigma
  ))

}

# EM algorithm implementation
lnm_EM = function(Y, Z, X, maxiter) {

  # get dimensions
  n = length(Y)
  q1 = ncol(Z)
  q2 = ncol(X)

  # generate initial parameter values
  # some random
  beta1_init = solve(t(Z) %*% Z) %*% t(Z) %*% Y # LS estimate
  beta2_init = rnorm(q1, 0, sd(Y))
  gamma_init = rnorm(q2, 0, 2) # values from roughly +-5
  sigma_init = sd(Y)

  # initial update
  result = update_EM(Y, Z, X, beta1_init, beta2_init, sigma_init, gamma_init)
  obj = lnm_logl(result$beta1, result$beta2, result$sigma, result$gamma, Y, Z, X)


  # main loop
  for (iter in 0:maxiter) {

    obj_old = obj

    # update
    result = update_EM(Y, Z, X, beta1, beta2, sigma, gamma)
    beta1 = result$beta1
    beta2 = result$beta2
    gamma = result$gamma
    sigma = result$sigma
    obj = lnm_logl(beta1, beta2, sigma, gamma, Y, Z, X)

    # print iteration number of objective value
    cat("Iter: ", iter, " obj: ", obj, "\n")

    # check convergence criterion
    ftolrel = 1e-12
    if (obj - obj_old < ftolrel*(abs(obj_old)+1))
      break


  }

  # return
  return(list(
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma,
    sigma = sigma
  ))
}


