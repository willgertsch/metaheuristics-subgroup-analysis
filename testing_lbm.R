# testing the logistic binomial mixture model

# create data set
set.seed(1234)
beta1 = c(-2, 0) # baseline rate of PrEP use: 0.1192029
beta2 = c(0, 1.2) # 0.3 rate in treatment and subgroup
gamma = c(-1.39, 1.79)
N = 1000
sex = rbinom(N, 1, 0.4) # female = 1
X = cbind(rep(1, N), sex)
treat = rep(c(0, 1), N/2)
Z = cbind(rep(1, N), treat)

# generate random trial sizes
sizes = sample(c(1,2,3,4,5,6), N, replace = T)
#sizes = rep(6, N)

# generate using function
mu1 = Z %*% (beta1 + beta2)
mu2 = Z %*% beta1
eta = X %*% gamma
gen_data = rlbm(N, mu1, mu2, eta, sizes)
Y = gen_data$Y
class = gen_data$class

df = data.frame(
  Y = Y,
  treat = as.factor(treat),
  sex = as.factor(sex),
  n = sizes,
  p = Y/sizes,
  class = class
)


# model without class knowledge
mod = glm(cbind(Y, n) ~ treat, family = binomial(), data = df)
summary(mod)

# model with class interactions
mod2 = glm(cbind(Y, n) ~ treat*class, family = binomial(), data = df)
summary(mod2)



# test likelihood
lbm_logl(beta1, beta2, gamma, Y, sizes, Z, X)

# fit model
library(metaheuristicOpt)
result = fit_lbm(Y, sizes, Z, X, 1000, 100, "PSO")
result


cbind(predict_class(X, result$gamma)[,1], class)
sum(predict_class(X, result$gamma)[,1] == class)/N

# calculated values
# estimated usage proportions
ilogit(t(c(1,0)) %*% result$beta1) # placebo
ilogit(t(c(1,1)) %*% result$beta1) # treatment
ilogit(t(c(1,0)) %*% (result$beta1+result$beta2)) # placebo, subgroup
ilogit(t(c(1,1)) %*% (result$beta1+result$beta2)) # treatment, subgroup

ilogit(t(c(1,0)) %*% beta1) # placebo
ilogit(t(c(1,1)) %*% beta1) # treatment
ilogit(t(c(1,0)) %*% (beta1+beta2)) # placebo, subgroup
ilogit(t(c(1,1)) %*% (beta1+beta2)) # treatment, subgroup

