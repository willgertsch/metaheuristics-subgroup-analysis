# simulation study
# compare performance of different algorithms
# key goal will be show that metaheuristics do better
# will compare class prediction performance and inference

# function for generating the LNM data
# returns list
generate_data = function(N, q2, beta1, beta2, sigma, p_class) {

  # generate q2 covariates
  # assuming scaled between 0 and 1
  X = matrix(runif(N * q2), ncol = q2)
  X = cbind(rep(1, N), X) # add intercept

  # randomly assign ~1/2 to treatment
  # should be fine to assign first half
  treat = c(rep(1, N/2), rep(0, N/2))
  Z = cbind(rep(1, N), treat)

  # randomly generate classes
  class = rbinom(N, 1, p_class)

  # get gamma from logistic regression
  mod = glm(class ~ X - 1, family = binomial())
  gamma = coef(mod)

  # compute means and response data
  mu1 = Z %*% (beta1 + beta2)
  mu2 = Z %*% beta1
  eta = X %*% gamma
  Y = rlnm(N, mu1, mu2, eta, sigma)


  sim_dat = list(
    Y = Y,
    X = X,
    Z = Z,
    class = class,
    gamma = gamma
  )

  return(sim_dat)

}

# test fitting model to this data
set.seed(1234)
beta1 = c(80, 0)
beta2 = c(0, 20)
sigma = 2
test = generate_data(1000, 10, beta1, beta2, sigma, 0.2)
result = fit_lnm(test$Y, test$Z, test$X, 1000, 100, "DE")
result

# EM
EM = lnm_EM(test$Y, test$Z, test$X, 100, silent = T)
EM

# simulation study
# main comparison metrics:
# bias of interaction effect estimates
# variance of interaction effect estimates
# rate of correct classification
# factors to vary:
# sample size (300, 500, 700, 1000)
# interaction effect size (10, 20, 30)
# latent factor rate (0, 0.2, 0.5)
# number of covariates => fix at 5
# algorithm

# think I just need to keep things simple and vary one factor at a time
# in order of importance
# algorithm
# latent factor rate
# interaction effect rate

################################################################################
# Simulation study 1
# how do different algorithms affect performance?
################################################################################
algorithms = c("DE", "MFO", "BHO", "HS")
algorithms = c(
  "PSO",
  #"ALO", slow
  "GWO",
 # "DA", slow
  "FFA",
  "GA",
  #"GOA", slow
  "HS",
  "MFO",
  "SCA",
  "WOA",
  "CLONALG",
  "DE",
  #"SFL", slow
  #"CSO", seems to broken
  #"ABC",
  #"KH", slow
  "CS",
  "BA",
  #"GBS",
  "BHO"
)
algorithms = c("GBS")

# data generating parameters
beta1 = c(80, 0)
beta2 = c(0, 20)
sigma = 2
N = 500
q2 = 10
p_rate = 0.2

# simulation parameters
samples = 10
iter = 1000
swarm = 100


results = rep(0, 7)
i = 1
for (s in 1:samples) {

  # generate data
  set.seed(i)
  test = generate_data(N, q2, beta1, beta2, sigma, p_rate)

  # fit all using all algorithms
  for (alg in algorithms) {

    cat("Sample ", s, "/", samples, " with ", alg, "\n")
    # fit
    mod = fit_lnm(test$Y, test$Z, test$X, iter, swarm, alg)

    # save
    p = sum(test$class == predict_class(test$X, mod$gamma)[,1])/N
    results = rbind(results,
                    c(
                      i,
                      mod$ll,
                      mod$beta1[1],
                      mod$beta1[2],
                      mod$beta2[1],
                      mod$beta2[2],
                      p
                      ))

  }
  i = i + 1
}

results = results[-1, ]
results = as.data.frame(results)
colnames(results) = c("seed", "ll", "beta11", "beta12", "beta21", "beta22", "p")
rownames(results) = NULL
results$algorithm = rep(algorithms, samples)

results %>%
  group_by(algorithm) %>%
  summarise(
    beta11_bias = abs(mean(beta11) - beta1[1]),
    beta11_se = sd(beta11),
    beta12_bias = abs(mean(beta12) - beta1[2]),
    beta12_se = sd(beta12),
    beta21_bias = abs(mean(beta21) - beta2[1]),
    beta21_se = sd(beta21),
    beta22_bias = abs(mean(beta22) - beta2[2]),
    beta22_se = sd(beta22),
    phat = mean(p)
  ) %>%
  arrange(desc(phat))


# repeat for EM algorithm ######################################################
EM_results = rep(0, 7)
i = 1
for (s in 1:samples) {

  # generate data
  set.seed(i)
  test = generate_data(N, q2, beta1, beta2, sigma, p_rate)

  # fit
  mod = lnm_EM(test$Y, test$Z, test$X, 100, silent = T)

  # save
  p = sum(test$class == predict_class(test$X, mod$gamma)[,1])/N
  EM_results = rbind(EM_results,
                  c(
                    i,
                    mod$ll,
                    mod$beta1[1],
                    mod$beta1[2],
                    mod$beta2[1],
                    mod$beta2[2],
                    p
                  ))
  i = i + 1
}

EM_results = EM_results[-1, ]
EM_results = as.data.frame(EM_results)
colnames(EM_results) = c("seed", "ll", "beta11", "beta12", "beta21", "beta22", "p")
rownames(EM_results) = NULL
EM_results$algorithm = rep("EM", samples)
