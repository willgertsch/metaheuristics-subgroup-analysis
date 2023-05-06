# testing for logistic normal model

################################################################################
# EM algorithm, M step for gamma, max_Q1 function
################################################################################
# like a logistic regression where only the proportions are known
# create test dataset
# from BMDS
df = data.frame(
  dose = c(0, 50, 100, 200, 400),
  N = c(20, 20, 20, 20, 20),
  incidence = c(0, 1, 2, 10, 19),
  p = c(0/20, 1/20, 2/20, 10/20, 19/20)
)
# example logistic model
mod = nls(p ~ 1/(1 + exp(-a-b*dose)), data = df, start = list(a=0, b=0))
plot(p ~ dose, data = df)
lines(predict(mod, type = "response") ~ df$dose)

# using NR algorithm
# beta = beta + (observed info matrix)^-1 * gradient

# data
W = df$p
X = cbind(rep(1, 5), df$dose)

# test function
gamma = max_Q1(X, W, c(0,0), 10)
gamma

# compare to model
coef(mod)
# my estimates are closer to what BMDS has

# plot
plot(p ~ dose, data = df)
lines(ilogit(X %*% gamma) ~ df$dose)


################################################################################
# EM algorithm, M step for betas, sigma, max_Q2 function
################################################################################
# a mixture regression where the proportions are already known

# create test dataset
# sample random weights
set.seed(1234)
n = 100
weights = rbeta(n, 1, 1)
z = runif(n)
Z = cbind(rep(1, n), z)
beta1 = c(80, 0)
beta2 = c(2, 20)
tbeta1 = beta1 + beta2
tbeta2 = beta1
sigma = 3

# generate from this distribution
Y = numeric(n)
for (i in 1:n) {

  class = rbinom(1, 1, weights[i])
  if (class == 1)
    Y[i] = rnorm(1, mean = Z[i, ] %*% tbeta1, sigma)
  else if (class == 0)
    Y[i] = rnorm(1, mean = Z[i,] %*% tbeta2, sigma)
}

# visualize
plot(Y ~ z)

# estimate the betas using
# shouldn't be using tbetas here as that is what I am trying to estimate
f1 = exp(-1/(2*sigma^2) * (Y - Z %*% tbeta1)^2)
f2 = exp(-1/(2*sigma^2) * (Y - Z %*% tbeta2)^2)
omega1 = 1/(weights * f1 + (1 - weights) * f2) * f1
omega2 = 1/(weights * f1 + (1 - weights) * f2) * f2
Omega1 = diag(c(omega1))
Omega2 = diag(c(omega2))

hat_tbeta1 = solve(t(Z) %*% Omega1 %*% Z) %*% t(Z) %*% Omega1 %*% Y
hat_tbeta2 = solve(t(Z) %*% Omega2 %*% Z) %*% t(Z) %*% Omega2 %*% Y
hat_tbeta2 # beta1
hat_tbeta1 - hat_tbeta2 # beta2
# huzzah it works

# estimate sigma
# use formula for estimating sigma_k then take average
sigma2_1 = sum(omega1*(Y - Z%*%hat_tbeta1)^2)/sum(omega1)
sigma2_2 = sum(omega2*(Y - Z%*%hat_tbeta2)^2)/sum(omega2)
hat_sigma = sqrt((sigma2_1 + sigma2_2)/2)
hat_sigma

# test function
beta1_init = c(mean(Y),0)
beta2_init = c(0,max(Y) - min(Y))
sigma_init = 1
result = max_Q2(Y, Z, weights, beta1_init, beta2_init, sigma_init)
result
# sequential update seems to work better

################################################################################
# Simulated data test case
################################################################################
# create dataset
set.seed(1234)
beta1 = c(80, 0)
beta2 = c(0, 30)
gamma = c(-1.39, 1.79)
sigma = 1
N = 100
sex = rbinom(N, 1, 0.4) # female = 1
X = cbind(rep(1, N), sex)
geneX = rbinom(100, 1, ilogit(X %*% gamma))
treat = rep(c(0, 1), N/2)
Z = cbind(rep(1, N), treat)

# generate using function
mu1 = Z %*% (beta1 + beta2)
mu2 = Z %*% beta1
eta = X %*% gamma
Y = rlnm(N, mu1, mu2, eta, sigma)

df = data.frame(
  Y = Y,
  treat = as.factor(treat),
  sex = as.factor(sex),
  geneX = as.factor(geneX)
)

# inspect
library(ggplot2)
ggplot(df, aes(x = Y, fill = treat)) +
  geom_histogram()
# looks like I did this correctly


# test likelihood
lnm_logl(beta1, beta2, sigma, gamma, Y, Z, X)

# fitting using metaheuristics
library(metaheuristicOpt)
result = fit_lnm(Y, Z, X, 1000, 100, "DE")

# test all algorithms in package
# takes a while to run
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

print(algorithms[1])
result_i = fit_lnm(Y, Z, X, 1000, 100, algorithms[1])
ll = lnm_logl(result_i$beta1, result_i$beta2, result_i$sigma, result_i$gamma,
              Y, Z, X)
results = c(unlist(result_i), ll, algorithms[1])
for (alg in algorithms[-1]) {
  print(alg)
  result_i = fit_lnm(Y, Z, X, 1000, 100, alg)
  ll = lnm_logl(result_i$beta1, result_i$beta2, result_i$sigma, result_i$gamma,
                Y, Z, X)
  results = rbind(results, c(unlist(result_i), ll, alg))
}
results
library(dplyr)
df = as.data.frame(results)
colnames(df)[8] = "ll"
colnames(df)[9] = "alg"
df %>% arrange(ll) %>% head(4)
# best: DE, MFO, BHO, HS

# EM algorithm
result = lnm_EM(Y, Z, X, 100)
result
# doesn't give as good result as for best metaheuristics

# prediction
# DE vs EM
class_DE = predict_class(X, c(-1.31431581238249, 1.86719412996483))
class_EM = predict_class(X, c(-1.306024, 1.858678))
# have the same counts
class_DE[,1] == class_EM[,1] # exact same prediction


# compare to geneX using classification rate
sum(class_DE[,1] == geneX)/n


# testing case where having singularity issues in EM
set.seed(7)
test = generate_data(500, 10, c(80,0), c(0,20), 2, 0.2)
mod = try(lnm_EM(test$Y, test$Z, test$X, 10, silent = T), FALSE)
length(mod)
if (length(mod) == 1) {
  print("The model has failed, write NA to data struct")
}

