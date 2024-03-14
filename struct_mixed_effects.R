# structured mixed effects model
# Shen and Qu (2020)


library(ggplot2)
library(metaheuristicOpt)
library(MASS)
library(tidyr)
library(dplyr)



################################################################################
# data generation
# multivariate normal strategy
################################################################################

# sample from MVN for each individual
set.seed(1234)

# covariance
generate_AR1_cov_matrix <- function(n, ar_coef) {
  # Create an empty covariance matrix
  cov_matrix <- matrix(0, nrow = n, ncol = n)

  # Populate the covariance matrix with AR(1) structure
  for (i in 1:n) {
    for (j in 1:n) {
      cov_matrix[i, j] <- ar_coef^abs(i - j)
    }
  }

  return(cov_matrix)
}

generate_data = function(
    N = 1000,
    timepoints = c(4,8,24,40),
    gene_prob = 0.3,
    alpha = c(-0.138, -0.086),
    mu1 = c(1.7, 2.2, 2.2, 2.1),
    mu0 = c(0.2, 0.4, 0.6, 0.7),
    gamma = c(-1, 2),
    sigma2 = 0.5,
    rho = 0.2
    ) {

  num_timepoints = length(timepoints)
  Z = cbind(rep(1, num_timepoints), timepoints)
  trt = rbinom(N, 1, 0.5)
  gene = rbinom(N, 1, gene_prob)

  D = sigma2 * diag(num_timepoints)
  R = generate_AR1_cov_matrix(num_timepoints, rho)

  result = matrix(data = NA, nrow = N, 2*num_timepoints+1)
  for (i in 1:N) {
    # compute subgroup membership
    eta_i= gamma[1] + gene[i]*gamma[2]
    p_i = exp(eta_i)/(1+exp(eta_i))
    delta_i = rbinom(1, 1, p_i)

    mu_i = mu1 * delta_i + mu0 * (1 - delta_i)
    big_mu_i = c(Z %*% alpha + trt[i]*mu_i, mu_i)

    Sigma_i = trt[i]^2 * D + sigma2*R
    big_cov_i = rbind(cbind(Sigma_i, trt[i]*D), cbind(trt[i]*D, D))

    result[i, ] = c(i, mvrnorm(
      n = 1,
      mu = big_mu_i,
      Sigma = big_cov_i
    ))

  }
  colnames(result) = c('ID','Y4', 'Y8', 'Y24', 'Y40', 'b4', 'b8', 'b24', 'b40')

  Y_time = as.data.frame(result) %>% dplyr::select(starts_with('Y'), 'ID') %>%
    pivot_longer(cols = starts_with('Y'),names_to = 'time', values_to = 'Y') %>%
    mutate(time = substring(time, 2))

  b_time = as.data.frame(result) %>% dplyr::select(starts_with('b'), 'ID') %>%
    pivot_longer(cols = starts_with('b'),names_to = 'time', values_to = 'b') %>%
    mutate(time = substring(time, 2))

  dat = inner_join(Y_time, b_time, by = c('ID', 'time')) %>%
    mutate(
      trt = rep(trt, each = 4),
      gene = rep(gene, each = 4),
      time = as.numeric(time)
    )

  return(dat)
}

dat = generate_data(N = 1000,
                    timepoints = c(4,8,24,40),
                    gene_prob = 0.3,
                    alpha = c(-0.138, -0.086),
                    mu1 = c(1.7, 2.2, 2.2, 2.1),
                    mu0 = c(0.2, 0.4, 0.6, 0.7),
                    gamma = c(-1, 2),
                    sigma2 = 0.5,
                    rho = 0.2
)

ggplot(dat, aes(x=time, y = Y, group = ID, color = as.factor(gene))) +
  geom_point() + geom_line() +
  facet_wrap(~trt)


# log likelihood function
ll_fn = function(
    theta
    ) {

  # decode parameters
  alpha = theta[1:2]
  mu1 = theta[3:6]
  mu0 = theta[7:10]
  gamma = theta[11:12]
  sigma2 = theta[13]
  rho = theta[14]

  R = generate_AR1_cov_matrix(length(mu1), rho)

  # compute log-likelihood by summing lpdfs
  ll = 0
  for (i in 1:N) {

    # data for subject i
    dat_i = dat %>% filter(ID == i)
    Y_i = dat_i$Y
    trt_i = dat_i$trt[1]
    gene_i = dat_i$gene[1]

    # covariance
    Sigma_i = sigma2*(diag(trt_i, length(mu1)) + R)
    Sigma_inv_i = solve(Sigma_i)

    # subgroup membership
    eta_i = gamma[1] + gene_i*gamma[2]
    w_i = exp(eta_i)/(1+exp(eta_i))

    ll = ll - 0.5*log(det(Sigma_i)) - 0.5*t(Y_i)%*%Sigma_inv_i%*%Y_i +
      log(w_i*exp(t(Y_i)%*%Sigma_inv_i%*%mu1 - 0.5*t(mu1)%*%Sigma_inv_i%*%mu1) +
            (1-w_i)*exp(t(Y_i)%*%Sigma_inv_i%*%mu0 - 0.5*t(mu0)%*%Sigma_inv_i%*%mu0))


  }

  ll = ll - 0.5*N*length(mu1)*log(2*pi)
  print(as.numeric(ll))
  return(as.numeric(ll))

}

N = 1000
dat
ll_fn(c(-0.138, -0.086, 1.7, 2.2, 2.2, 2.1, 0.2, 0.4, 0.6, 0.7, -1, 2, 0.5, 0.2))

lower = c(-1, -1, # alpha
          0, 0, 0, 0, # mu1
          0, 0, 0, 0, # mu0
          -2, 0, # gamma
          0.1, 0.1 # sigma2, rho
          )
upper = c(1, 1, # alpha
          5, 5, 5, 5, # mu1
          3, 3, 3, 3, # mu0
          0, 5, # gamma
          1, 1 # sigma2, rho
)

numVar = 14
rangeVar = rbind(lower,upper)

out = metaOpt(ll_fn, optimType = "MAX", algorithm = "HS", numVar=14, rangeVar=rangeVar,
        control = list(), seed = NULL)
