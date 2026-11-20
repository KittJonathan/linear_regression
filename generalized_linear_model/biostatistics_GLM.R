# Generalized Linear Model
# Youtube playlist: https://youtube.com/playlist?list=PLl4b1XM725Oc8g3V9p98UeRTOXoAlGEYw&si=y8mYeSswR-6WgYJW
# Github repo: https://github.com/FransRodenburg/Biostatistics/tree/main/GLM

# Packages ----

library(tidyverse)

theme_set(theme_bw())

# Data ----

apples <- read_csv("generalized_linear_model/example_data_counts.csv")

# Count regression ----

apples |> 
  ggplot(aes(x = age, y = apples)) +
  geom_point()

# Poisson distribution: how does lambda (average nb of apples) depend on x (age)?

# Simple linear model:
m1 <- lm(apples ~ age, data = apples)
m1$coefficients

apples |> 
  ggplot(aes(x = age, y = apples)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

plot(m1, which = 1)
plot(m1, which = 2)
plot(m1, which = 3)
plot(m1, which = 5)

## GLM ----

GLM <- glm(apples ~ age, family = poisson(link = "log"), data = apples)

# Enter the estimates
intercept   <- coef(GLM)[1]
slope       <- coef(GLM)[2]
inverselink <- exp

# Create a large sequence of numbers along the x-axis
x <- seq(0, 3, length.out = 1000)

# Compute the predictions on the scale of the linear predictor
eta <- intercept + slope * x

# Compute the resulting values on the y-axis
y <- inverselink(eta)

# Plot the resulting model
plot(x, y, type = "l")

# The following gives us two things:
#  - Predictions of lambda for every value of x
#  - A standard error at that value of x
pred <- predict(GLM, newdata = list(age = x), type = "link", se.fit = TRUE)

# Combined, those two can create a 95% CI on the scale of linear predictor.
# On this scale (not on the original scale!), we can use a normal-approximation.
eta     <- pred$fit
eta_lwr <- pred$fit + qnorm(0.025) * pred$se.fit
eta_upr <- pred$fit + qnorm(0.975) * pred$se.fit

# Compute the resulting values on the y-axis
y     <- inverselink(eta)
y_lwr <- inverselink(eta_lwr)
y_upr <- inverselink(eta_upr)

# Plot the resulting model
plot(x, y, type = "l")

# Add the confidence band
lines(x, y_lwr, lty = "dashed")
lines(x, y_upr, lty = "dashed")

# The following gives us two things:
#  - Predictions of lambda for every value of x
#  - A standard error at that value of x
pred <- predict(GLM, newdata = list(age = x), type = "link", se.fit = TRUE)

# Combined, those two can create a 95% CI on the scale of linear predictor.
# On this scale (not on the original scale!), we can use a normal-approximation.
eta     <- pred$fit
eta_lwr <- pred$fit + qnorm(0.025) * pred$se.fit
eta_upr <- pred$fit + qnorm(0.975) * pred$se.fit

# Compute the resulting values on the y-axis
# For the prediction interval, we want the lowest 2.5% of a Poisson distribution
# with lambda equal to the lower bound of the confidence interval, and the 
# highest 2.5% o a Poisson distribution with lambda equal to the upper bound.
y     <- inverselink(eta)
y_lwr <- qpois(0.025, lambda = inverselink(eta_lwr))
y_upr <- qpois(0.975, lambda = inverselink(eta_upr))

# Plot the resulting model
plot(x, y, type = "l")

# Add the prediction band
lines(x, y_lwr, lty = "dashed")
lines(x, y_upr, lty = "dashed")
