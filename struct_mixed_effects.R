# structured mixed effects model
# Shen and Qu (2020)

library(dplyr)
library(ggplot2)
library(metaheuristicOpt)
library(MASS)
library(tidyr)

# data generation
# SEM strategy
# assuming balanced data, so no unique covariance matrix
set.seed(1234)
N = 1000
time = rep(c(4,8,24,40), N)
timepoints = length(time)/N
ID = rep(1:N, each = timepoints)
trt = rep(rbinom(N, 1, 0.5), each=timepoints)
Z = cbind(rep(1, N*timepoints), time)

# subgroup predictor
gene = rbinom(N, 1, 0.2)
X = cbind(rep(1, N), gene)
gamma = c(-1, 2)
p = exp(X %*% gamma)/(1 + exp(X %*% gamma))

# generate from mixture
# have to transpose and vec
mu1 = c(1.7, 2.2, 2.2, 2.1)
mu0 = c(0.2, 0.4, 0.6, 0.7)
D = 1*diag(timepoints)
subgroup = rbinom(N, 1, p)
b = mvrnorm(N, mu0, D)
for (i in nrow(b)) {
  if (subgroup[i])
    b[i, ] = mvrnorm(1, mu1, D)
}
epsilon = c(t(mvrnorm(N, rep(0, timepoints), 0.428*diag(timepoints))))
alpha = c(-0.138, -0.086)
Y = Z %*% alpha + trt * c(t(b)) + epsilon

dat = data.frame(
  ID,
  trt,
  time,
  Y,
  gene = rep(gene, each=timepoints)
)
dat

ggplot(dat, aes(x=time, y=Y, group=ID, color=as.factor(trt))) +
  geom_point() + geom_line() +
  facet_wrap(~gene)

# data generation
# multivariate normal strategy

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
    delta_i = rbinom(1, 1, p)

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

ggplot(dat, aes(x=time, y = Y, group = ID, color = as.factor(trt))) +
  geom_point() + geom_line() +
  facet_wrap(~gene)


# log likelihood function
ll_fn = function(
    theta
    ) {

  # decode parameters
  alpha = theta[1:2]
  mu1 = theta[3:7]
  mu0 = theta[8:11]
  gamma = theta[12:13]
  sigma2 = theta[14]
  rho = theta[15]

  # compute log-likelihood by summing lpdfs
  ll = 0
  for (i in 1:N) {

    Y_i = dat %>% filter(ID == i) %>%
      .$Y

    big_mu_i = # problem: its only normal conditional on unobserved delta_i
      # need to use law of total probability

    ll = ll -0.5*log(det(big_sigma_i)) -
      0.5 * t(Y_i - big_mu_i) %*% inv(big_sigma_i) %*% (Y_i - big_mu_i)
  }

}


ll_fn(c(-0.138, -0.086, 1.7, 2.2, 2.2, 2.1, 0.2, 0.4, 0.6, 0.7, -1, 2, 0.5, 0.2))
