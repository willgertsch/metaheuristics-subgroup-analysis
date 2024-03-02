# experiment with using metaheuristics to identify subgroups
# the idea is use nature-inspired algorithms instead of a tree
# Nature Inspired Subgroup Search

# libraries
library(lme4)
library(ecr)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(shiny)

# generate longitudinal data with a subgroup effect
# N: sample size
# timepoints:
# b0: fixed baseline mean
# b1: treatment difference at baseline
# b2: time effect
# b3: overall trt:time effect
# b4: trt time effect modified by subgroup
# b_se: standard error of random intercept
generate_data = function(N = 1000,
                         timepoints = 6,
                         outcome = 'normal',
                         b0 = 0.5,
                         b_se = 0.3,
                         b1AB = 0,
                         b1AC = 0,
                         b1ABC = 0,
                         b2 = 0.5,
                         b3AB = 0.2,
                         b3AC = 0.3,
                         b3ABC = 0.6,
                         b4AB = 0,
                         b4AC = 0,
                         b4ABC = 0
                         ) {

  ID = rep(1:N, each=timepoints)
  time = rep(0:(timepoints-1), N)
  #trt = rep(rbinom(N, 1, 0.5), each = timepoints)
  trt = rep(sample(c('A', 'A+B', 'A+C', 'A+B+C'), size = N, prob = rep(0.25, 4),
                   replace=T), each = timepoints) %>%
    factor(levels = c('A', 'A+B', 'A+C', 'A+B+C'))
  trtAB = trt == 'A+B'
  trtAC = trt == 'A+C'
  trtABC = trt == 'A+B+C'

  # generate baseline covariates
  X = data.frame(
    ID = 1:N,
    Z = rnorm(N,0, 1),
    N1 = rnorm(N, 50, 20),
    T1 = rt(N, 6),
    E1 = rexp(N, 1),
    B1 = rbinom(N, 1, 0.2),
    B2 = rbinom(N, 1, 0.4),
    B3 = rbinom(N, 1, 0.6)
  )

  # subgroup definition
  S = rep(ifelse(X$B1 == 1, 1, 0), each= timepoints)

  b = rep(rnorm(N, 0, b_se), each = timepoints) # random intercept

  # mean function with random intercept
  #mu =  b0 + b + (b1)*trt + (b2)*time + (b3+b4*S)*time*trt
  mu = b0 + b +
    b1AB*trtAB + b1AC*trtAC + b1ABC*trtABC +
    b2*time +
    (b3AB+b4AB*S)*trtAB*time + (b3AC+b4AC*S)*trtAC*time + (b3ABC+b4ABC*S)*trtABC*time

  # outcome
  if (outcome == 'normal')
    Y = rnorm(N*timepoints, mean = mu, sd = 1.5)
  else if (outcome == 'binary')
    Y = rbinom(N*timepoints, size = 1, prob = exp(mu)/(1+exp(mu)))
  else if (outcome == 'count')
    Y = rpois(N*timepoints, exp(mu))

  d = data.frame(
    ID=ID,
    Y=Y,
    time=time,
    trt=trt,
    S=S
  ) %>%
    left_join(X, by = 'ID')

  return(d)
}

# generate longitudinal dataset
set.seed(12340)
# this estimates are based on PrEP model
d = generate_data(N = 1000,
                  timepoints = 6,
                  outcome = 'binary',
                  b0 = -4.8,
                  b_se = 3.6,
                  b1AB = 0.5,
                  b1AC = 0.4,
                  b1ABC = 0.43,
                  b2 = 0.018,
                  b3AB = -0.046,
                  b3AC = 0.0009,
                  b3ABC = 0.19,
                  b4AB = 0,
                  b4AC = 0,
                  b4ABC = 0
)

# plot by treatment arm
# individual effects
ggplot(data = d, aes(x=time, y=Y, group=ID,color=trt)) +
  geom_line()

# population effects
d %>%
  group_by(trt, time) %>%
  summarise(mu = mean(Y)) %>%
  ggplot(aes(x=time, y=mu, color=trt))+
  geom_point()+geom_line()

# subgroup
# individual effects
ggplot(data = d, aes(x=time, y=Y, group=ID,color=as.factor(trt))) +
  geom_line() +
  facet_wrap(~S)

# population effects
d %>%
  group_by(trt, time, S) %>%
  summarise(mu = mean(Y)) %>%
  ggplot(aes(x=time, y=mu, color=trt))+
  geom_point()+geom_line() +
  facet_wrap(~S)

# modeling
mod = glmer(Y ~ (1|ID) + trt*time, data = d, family = binomial, nAGQ = 25)
summary(mod)

# include subgroup interactions
mod_S = glmer(Y ~ (1|ID) + trt*time*S, data = d, family = binomial, nAGQ = 25)
ss = summary(mod_S)
ss$coefficients

anova(mod, mod_S)


# shiny app to experiment with
# Define UI for application
ui <- fluidPage(
  titlePanel("Generate Data and Plot"),

  sidebarLayout(
    sidebarPanel(
      sliderInput("N", "Number of data points:", value = 1000, min = 10, max = 5000),
      sliderInput("timepoints", "Number of timepoints:", value = 6, min = 2, max = 20),
      sliderInput("b0", "b0:", value = 0.5, min = -10, max = 10, step = 0.1),
      sliderInput("b1", "b1:", value = 0, min = -10, max = 10, step = 0.1),
      sliderInput("b2", "b2:", value = 0.5, min = -10, max = 10, step = 0.1),
      sliderInput("b3", "b3:", value = 1, min = -10, max = 10, step = 0.1),
      sliderInput("b4", "b4:", value = 0.5, min = -10, max = 10, step = 0.1),
      sliderInput("b_se", "b_se:", value = 0.3, min = .0001, max = 10, step = 0.1)
    ),

    mainPanel(
      plotOutput("dataPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {

  # Reactive expression to generate data
  d <- reactive({
    generate_data(N = input$N, timepoints = input$timepoints,
                  b0 = input$b0, b1 = input$b1, b2 = input$b2,
                  b3 = input$b3, b4 = input$b4, b_se = input$b_se)
  })

  # Render the plot
  output$dataPlot <- renderPlot({
    ggplot(data = d(), aes(x=time, y=Y, group=ID,color=as.factor(trt))) +
      geom_line() +
      facet_wrap(~S)
  })
}

# Run the application
shinyApp(ui = ui, server = server)


# find subgroup using nature inspired algorithm
# define fitness function
# gene is a binary vector that says a binary covariate should be including in
# subgroup definition
# B1 AND B2 AND B3 AND .....
# gene controls inclusion in this AND statement
X = d %>%
  distinct(ID, .keep_all = T) %>%
  select(starts_with('B')) %>%
  as.matrix()
fn = function(gene) {

  gene = as.logical(gene)
  if (sum(gene) == 0)
    return(-Inf)

  if (is.vector(X[, gene]))
    S1 = X[, gene] == 1
  else if (is.matrix(X[, gene]))
    S1 = rowSums(X[, gene]) == sum(gene)

  S1 = rep(S1, each = 6)
  mod = lmer(Y ~ (1|ID) + trt*time*S1, data = d)
  ss = summary(mod)
  abs(ss$coefficients[8,4])
}

fn(c(1, 1, 0))

res = ecr(
  fn,
  minimize = F,
  n.bits = 3L,
  n.objectives = 1L,
  representation = "binary",
  mu = 10L, # population size
  lambda = 5, # number of individuals generated in each generation,
)
res
res$best.x
res$last.population
res$message
