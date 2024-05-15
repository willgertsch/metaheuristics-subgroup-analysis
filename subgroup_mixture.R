# generalized subgroup mixtures

# generate data
# N: sample size
# beta0: main effects
# beta1: subgroup interaction effects
# gamma: subgroup membership covariates
# type: outcome type, "continuous", "binary", 'count"
# sigma: for continuous outcome, set the two variance parameters
# trt_prob: treatment assignment probability
# X: subgroup covariate data.frame turned into matrix by model.matrix
generate_data = function(N, beta0, beta1, gamma, type, sigma = NULL, trt_prob, X) {


  # sanitity checks
  if (length(beta0) != length(beta1))
    stop('beta0, beta1 dimension mismatch')


  # design matrix
  Z = cbind(rep(1, N), rbinom(N, 1, trt_prob))

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

predict_class = function(X, gamma) {

  eta = X %*% gamma
  p = ilogit(eta)
  class = ifelse(p > 0.5, 1, 0)
  return(cbind(class, p))

}
