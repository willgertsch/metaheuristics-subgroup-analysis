# simulation study
# compare performance of different algorithms
# key goal will be show that metaheuristics do better
# will compare class prediction performance and inference

# function for generating the LNM data
# returns list
generate_data = function(N, p, q2) {

  beta1 = c(rnorm(1, 80, 15), rnorm(1, 0, 1)) # assuming no intervention effect
  beta2 = c(rnorm(1, 0, 1), rnorm(1, 30, 10)) # intervention effect in subgroup
  sigma = runif(1, 0.5, 10)

  # generate classes
  class = rbinom(N, 1, p)

  # generate q2 covariates
  # assuming scaled between 0 and 1
  X = matrix(runif(N * q2), ncol = q2)
  X = cbind(rep(1, N), X) # intercept

  # fit a logistic regression model to get gamma estimates
  mod = glm(class ~ X - 1, family = binomial())
  gamma = coef(mod)

  # randomly assign ~1/2 to treatment
  treat = rbinom(N, 1, 0.5)
  Z = cbind(rep(1, N), treat)

  mu1 = Z %*% (beta1 + beta2)
  mu2 = Z %*% beta1
  eta = X %*% gamma
  Y = rlnm(N, mu1, mu2, eta, sigma)

  sim_dat = list(
    Y = Y,
    X = X,
    Z = Z,
    class = class,
    beta1 = beta1,
    beta2 = beta2,
    gamma = gamma,
    sigma = sigma
  )

  return(sim_dat)

}

# test fitting model to this data
set.seed(1234)
test = generate_data(1000, 0.2, 10)
result = fit_lnm(test$Y, test$Z, test$X, 1000, 100, "DE")
result

# EM
EM = lnm_EM(test$Y, test$Z, test$X, 100, silent = T)
EM

# compare against true values using L2 norm
sum(result$beta1 - test$beta1)^2
sum(EM$beta1 - test$beta1)^2

sum(result$beta2 - test$beta2)^2
sum(EM$beta2 - test$beta2)^2

sum(result$gamma - test$gamma)^2
sum(EM$gamma - test$gamma)^2

(result$sigma - test$sigma)^2
(EM$sigma - test$sigma)^2

sum(predict_class(test$X, result$gamma)[,1] == test$class)/1000
sum(predict_class(test$X, EM$gamma)[,1] == test$class)/1000

################################################################################
# Simulation study 1
# how do different algorithms affect performance?
################################################################################
#algorithms = c("DE", "MFO", "BHO", "HS")
algorithms = c(
  "PSO",
  "ALO",
  "GWO",
  "DA",
  "FFA",
  "GA",
  "GOA",
  "HS",
  "MFO",
  "SCA",
  "WOA",
  "CLONALG",
  "DE",
  "SFL",
  #"CSO", seems to broken
  #"ABC",
  "KH",
  "CS",
  "BA",
  #"GBS",
  "BHO"
)
samples = 10
N = 1000
p = 0.2
q = 10
set.seed(4321)

results = numeric(6)

for (s in 1:samples) {

  # generate data
  test = generate_data(N, p, q)

  for (alg in algorithms) {
    cat("Sample ", s, " of ", samples, " using ", alg, "\n")

    # fit model using algorithm
    out = fit_lnm(test$Y, test$Z, test$X, 1000, 100, alg)

    # save
    results_i = c(
      lnm_logl(out$beta1, out$beta2, out$sigma, out$gamma,
               test$Y, test$Z, test$X),
      sum(out$beta1 - test$beta1)^2,
      sum(out$beta2 - test$beta2)^2,
      sum(out$gamma - test$gamma)^2,
      (out$sigma - test$sigma)^2,
      sum(predict_class(test$X, out$gamma)[,1] == test$class)/N
    )

    results = rbind(results, results_i)
  }


  # EM algorithm
  out = lnm_EM(test$Y, test$Z, test$X, 100, silent = T)
  results_i = c(
    lnm_logl(out$beta1, out$beta2, out$sigma, out$gamma,
             test$Y, test$Z, test$X),
    sum(out$beta1 - test$beta1)^2,
    sum(out$beta2 - test$beta2)^2,
    sum(out$gamma - test$gamma)^2,
    (out$sigma - test$sigma)^2,
    sum(predict_class(test$X, out$gamma)[,1] == test$class)/N
  )

  results = rbind(results, results_i)

}

# process results
results = results[-1, ]
results = as.data.frame(results)
colnames(results) = c('ll', 'beta1', 'beta2', 'gamma', 'sigma', 'accuracy')
rownames(results) = NULL
results$algorithm = rep(c('EM', algorithms), samples)

# summarise
library(dplyr)
results %>%
  group_by(algorithm) %>%
  summarise(
    ll = median(ll),
    beta1 = median(beta1),
    beta2 = median(beta2),
    gamma = median(gamma),
    sigma = median(sigma),
    accuracy = median(accuracy)
  ) %>%
  arrange(desc(ll))
