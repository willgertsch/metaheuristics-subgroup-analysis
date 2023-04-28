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

# generate using function
mu1 = Z %*% (beta1 + beta2)
mu2 = Z %*% beta1
eta = X %*% gamma
geneX = ifelse(ilogit(eta) >= 0.5, 1, 0)
Y = rlbm(N, mu1, mu2, eta, sizes)

df = data.frame(
  Y = Y,
  treat = as.factor(treat),
  sex = as.factor(sex),
  n = sizes,
  p = Y/sizes,
  geneX = geneX
)

df %>%
  group_by(treat) %>%
  summarise(phat = mean(p))

mod = glm(cbind(Y, n) ~ treat, family = binomial(), data = df)
summary(mod)
# no significant treatment effect
# some times is significant

# stratify by geneX
summary(glm(cbind(Y, n) ~ treat, family = binomial(), data = df, subset = geneX==1))


# test likelihood
lbm_logl(beta1, beta2, gamma, Y, sizes, Z, X)

# fit model
library(metaheuristicOpt)
result = fit_lbm(Y, sizes, Z, X, 1000, 100, "HS")
result

# calculated values
# estimated usage proportions
ilogit(t(c(1,0)) %*% result$beta1) # placebo
ilogit(t(c(1,1)) %*% result$beta1) # treatment
ilogit(t(c(1,0)) %*% (result$beta1+result$beta2)) # treatment, subgroup
ilogit(t(c(1,1)) %*% (result$beta1+result$beta2)) # treatment, subgroup

ilogit(t(c(1,0)) %*% beta1) # placebo
ilogit(t(c(1,1)) %*% beta1) # treatment
ilogit(t(c(1,0)) %*% (beta1+beta2)) # treatment, subgroup
ilogit(t(c(1,1)) %*% (beta1+beta2)) # treatment, subgroup
