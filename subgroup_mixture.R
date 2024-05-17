# generalized subgroup mixtures

# generate data
# N: sample size
# beta0: main effects
# beta1: subgroup interaction effects
# beta can have 2 or 3 parameters
# gamma: subgroup membership covariates
# type: outcome type, "continuous", "binary", 'count"
# sigma: for continuous outcome, set the two variance parameters
# trt_prob: treatment assignment probability
# X: subgroup covariate data.frame turned into matrix by model.matrix
generate_data = function(N, beta0, beta1, gamma, type, sigma = NULL, trt_prob, X) {


  # sanitity checks
  if (length(beta0) != length(beta1))
    stop('beta0, beta1 dimension mismatch')


  baseline_var = rnorm(N)
  # design matrix
  # match Z to dimension of beta
  if (length(beta0) == 2) {
    Z = cbind(rep(1, N), rbinom(N, 1, trt_prob))
  }
  else if (length(beta0) == 3) {
    Z = cbind(rep(1, N), rbinom(N, 1, trt_prob), baseline_var)
    X = cbind(X, baseline_var)
  }
  else
    stop("generate_data: length(beta) > 3 adding more covariates is better in real life, but for simplicity we don't.")


  # compute means
  mu0 = Z %*% beta0
  mu1 = Z %*% (beta0 + beta1)


  eta = X %*% gamma
  class = rbinom(N, 1, ilogit(eta))

  Y = numeric(N)
  if (type == 'continuous') {
    for (i in 1:N) {
      if (class[i] == 0) {
        Y[i] = rnorm(1, mu0[i], sigma[1])
      }
      else if (class[i] == 1) {
        Y[i] = rnorm(1, mu1[i], sigma[2])
      }
    }
  }
  else if (type == 'binary') {
    for (i in 1:N) {
      if (class[i] == 0) {
        Y[i] = rbinom(1, 1, ilogit(mu0[i]))
      }
      else if (class[i] == 1) {
        Y[i] = rbinom(1, 1, ilogit(mu1[i]))
      }
    }
  }
  else if (type == 'count') {

    for (i in 1:N) {
      if (class[i] == 0) {
        Y[i] = rpois(1, exp(mu0[i]))
      }
      else if (class[i] == 1) {
        Y[i] = rpois(1, exp(mu1[i]))
      }
    }
  }

  return(
    list(
      Y=Y, Z=Z, X=X, class=class
    )
  )

}

# generate subgroup covariate data frame
# converts automatically into matrix form
# N: sample size
generate_X = function(N) {
  X1 = rbinom(N, 1, 0.5)
  X2 = rnorm(N)

  X = data.frame(X1, X2)

  temp = rep(1, nrow(X))
  X = model.matrix(temp ~., data = cbind(temp, X))

  return(X)
}


ilogit = function(eta) {
  1/(1 + exp(-eta))
}

# generate log-likelihood function for use with optimizer
#
ll_factory = function(type, Y, X, Z) {

  if (type == 'continuous') {
    ll_fun = function(var) {
      beta0 = var[1:2]
      beta1 = var[3:4]
      gamma = var[5:(4+ncol(X))]
      sigma = var[(5+ncol(X)):length(var)]

      # check constraint
      if (beta1[2] <= 0)
        return(-Inf)

      eta = X %*% gamma
      p = ilogit(eta)
      mu0 = Z %*% beta0
      mu1 = Z %*% (beta0 + beta1)
      f0 = dnorm(Y, mu0, sigma[1])
      f1 = dnorm(Y, mu1, sigma[2])
      ll = sum(log(p * f1 + (1 - p) * f0))
      if (is.na(ll))
        return(-Inf)
      else
        return(ll)
    }
  }
  else if (type == 'binary') {

    ll_fun = function(var) {

      beta0 = var[1:2]
      beta1 = var[3:4]
      gamma = var[5:(4+ncol(X))]

      # check constraint
      if (beta1[2] <= 0)
        return(-Inf)

      eta = X %*% gamma
      p = ilogit(eta)
      mu0 = Z %*% beta0
      mu1 = Z %*% (beta0 + beta1)
      f0 = dbinom(Y, 1, ilogit(mu0))
      f1 = dbinom(Y, 1, ilogit(mu1))
      ll = sum(log(p * f1 + (1 - p) * f0))
      if (is.na(ll))
        return(-Inf)
      else
        return(ll)
    }
  }
  else if (type == 'count') {
    ll_fun = function(var) {

      beta0 = var[1:2]
      beta1 = var[3:4]
      gamma = var[5:(4+ncol(X))]

      # check constraint
      if (beta1[2] <= 0)
        return(-Inf)

      eta = X %*% gamma
      p = ilogit(eta)
      mu0 = Z %*% beta0
      mu1 = Z %*% (beta0 + beta1)
      f0 = dpois(Y, exp(mu0))
      f1 = dpois(Y, exp(mu1))
      ll = sum(log(p * f1 + (1 - p) * f0))
      if (is.na(ll))
        return(-Inf)
      else
        return(ll)
    }
  }
  else
    stop('ll_factory: type not supported')

  return(ll_fun)

}

# model fitting function
# type: continuous, binary, count
# Y: outcome
# X: subgroup design matrix
# Z: outcome model design matrix
# swarm: population size for algorithm
# maxIter: maximum number of iterations
# algorithm: algorithm to use with metaOpt
# rangeVar: parameter bounds, can manually set otherwise default bounds will be set
fit_model = function(type, Y, X, Z, swarm = 40, maxIter = 500, algorithm, rangeVar=NULL) {

  # get likelihood function
  ll_fun = ll_factory(type, Y, X, Z)

  # set variable bounds
  # setting reasonable bounds for logistic variables.
  # basing this on Bayesian priors for logistic regression
  # http://www.stat.columbia.edu/~gelman/research/published/priors11.pdf
  logOR_lb = -10
  logOR_ub = 5 # more than a 0.01 to 0.5 jump, which is not going to happen for this data
  if (is.null(rangeVar)) {
    if (type == 'continuous') {
      minY = min(Y)
      maxY = max(Y)
      rangeY = maxY - minY
      sdY = sd(Y)
      rangeVar = matrix(c(
      minY, -rangeY, -rangeY, 0, rep(logOR_lb, ncol(X)), 0.1, 0.1,
      maxY, rangeY, rangeY, rangeY, rep(logOR_ub, ncol(X)), sdY, sdY
      ),
      nrow = 2, byrow = T)
    }
    else if (type == 'binary') {
      rangeVar = matrix(c(
        logOR_lb, logOR_lb, logOR_lb, 0, rep(logOR_lb, ncol(X)),
        logOR_ub, logOR_ub, logOR_ub, logOR_ub, rep(logOR_ub, ncol(X))
      ),
      nrow = 2, byrow = T)
    }
    else if (type == 'count') {

      # think the OR bounds will also be find for Poisson regression
      rangeVar = matrix(c(
        logOR_lb, logOR_lb, logOR_lb, 0, rep(logOR_lb, ncol(X)),
        logOR_ub, logOR_ub, logOR_ub, logOR_ub, rep(logOR_ub, ncol(X))
      ),
      nrow = 2, byrow = T)
    }
  }


  result = metaheuristicOpt::metaOpt(
    ll_fun,
    optimType = 'MAX',
    algorithm = algorithm,
    numVar = ncol(rangeVar),
    rangeVar = rangeVar,
    control = list(
      numPopulation = swarm,
      maxIter = maxIter
    )
  )

  # extract results
  param = as.numeric(result$result)
  beta1 = param[1:2]
  beta2 = param[3:4]
  gamma = param[5:(5 + ncol(X) - 1)]
  if (type == 'continuous')
    sigma = tail(param, 2)
  else
    sigma = NULL


  return(list(
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma,
    sigma = sigma,
    algorithm = algorithm,
    swarm = swarm,
    maxIter = maxIter,
    ll = as.numeric(result$optimumValue),
    timeElapsed = result$timeElapsed,
    Y = Y,
    X = X,
    Z = Z,
    type = type
  ))

}

# predict subgroup class based on model fit
predict_class = function(mod) {

  eta = mod$X %*% mod$gamma
  p = ilogit(eta)
  class = ifelse(p > 0.5, 1, 0)
  return(data.frame(class, p))

}

# nicely summarizes result from model
# mod: list output from fit_model
summary.subgroup = function(mod) {

  cat(mod$type, "logistic-mixture model fit using", mod$algorithm, '\n')
  cat("Log-likelihood:", mod$ll, '\n')
  cat("Estimated treatment effect for s=0:", mod$beta1[2], '\n')
  cat("Estimated treatment effect for s=1:", mod$beta1[2] + mod$beta2[2], '\n')
  cat('Treatment effect increase:', mod$beta2[2], '\n')

  class_dat = predict_class(mod)
  cat('Subgroup:', table(class_dat$class), '\n')

}
