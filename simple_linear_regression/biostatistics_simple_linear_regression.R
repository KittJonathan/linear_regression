# Simple Linear Regression
# Youtube playlist: https://www.youtube.com/playlist?list=PLl4b1XM725OfDg-0FeWlcF2sBTAveFBwY
# Github repo: https://github.com/FransRodenburg/Biostatistics/tree/main/SimpleLinearRegression

# Packages ----

library(tidyverse)
library(tidymodels)
library(ggpubr)

theme_set(theme_bw())

# Import data ----

df <- read_csv("simple_linear_regression/example_data_simple_linear_regression.csv")

head(df)

# Plot data ----

df |> 
  ggplot(aes(x = weight, y = height)) +
  geom_point()

# Regression line ----

linear_model <- lm(height ~ weight, data = df)

linear_model |> 
  augment() |> 
  ggplot(aes(x = weight, y = height)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE,
              color = "steelblue2", linewidth = 0.5) +
  geom_segment(aes(x = weight, xend = weight,
                   y = height, yend = .fitted),
               color = "steelblue", linewidth = 0.4, linetype = "dashed") +
  stat_regline_equation(label.x = 100, label.y = 155, decreasing = F) +
  coord_fixed()

linear_model |> 
  augment() |> 
  mutate(x1 = weight - .resid/2, 
         x2 = weight + .resid/2,
         y1 = height, 
         y2 = .fitted) |> 
  ggplot(aes(x = weight, y = height)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE,
              color = "steelblue2", linewidth = 0.5) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
            fill = "red", alpha = 0.1, color = "red", linewidth = 0.2) +
  stat_regline_equation(label.x = 100, label.y = 155, decreasing = F) +
  coord_fixed()

# Regression table ----

linear_model

summary(linear_model)

summary(linear_model)$coefficients

ggplot() +
  geom_point(data = df,
             aes(x = weight, y = height)) +
  geom_abline(intercept = linear_model$coefficients[[1]],
              slope = linear_model$coefficients[[2]]) +
  coord_fixed()

# Residuals ----

linear_model$residuals
summary(linear_model$residuals)

# Coefficients ----

linear_model$coefficients
summary(linear_model)

# Residual standard error
rse <- sqrt(sum(linear_model$residuals^2) / (length(linear_model$residuals) - 2))
all.equal(rse, summary(linear_model)$sigma)

# Standard error for the slope
rse / sqrt(sum((df$weight - mean(df$weight))^2))

# Standard error for the intercept
rse * sqrt((1 / length(df$weight)) + ((mean(df$weight))^2 / sum((df$weight - mean(df$weight))^2)))

x <- df$weight
y <- df$height
n <- length(x)

# Residual standard error
s <- sqrt(sum(resid(linear_model)^2) / (n - 2))

Sxx <- sum((x - mean(x))^2)

se_slope     <- s / sqrt(Sxx)
se_intercept <- s * sqrt(1/n + mean(x)^2 / Sxx)

c(se_intercept = se_intercept, se_slope = se_slope)

X <- model.matrix(linear_model)          # design matrix (column of 1s + x)
s2 <- sum(resid(linear_model)^2) / df.residual(linear_model)
vcov_manual <- s2 * solve(t(X) %*% X)
sqrt(diag(vcov_manual))

sqrt(diag(vcov(linear_model)))
coef(summary(linear_model))[, "Std. Error"]

est <- coef(linear_model)
se  <- c(se_intercept, se_slope)
t   <- est / se
p   <- 2 * pt(-abs(t), df = n - 2)
cbind(est, se, t, p)

# t-value = (estimate - 0) / standard error

(summary(linear_model)$coefficients[1] - 0) / summary(linear_model)$coefficients[3]
summary(linear_model)$coefficients[5]

(summary(linear_model)$coefficients[2] - 0) / summary(linear_model)$coefficients[4]
summary(linear_model)$coefficients[6]

# pt(q = summary(linear_model)$coefficients[5], df = 29, lower.tail = F)


# Setup

X <- model.matrix(linear_model)               # design matrix (includes intercept column)
y <- model.response(model.frame(linear_model))
n <- nrow(X)                                  # number of observations
p <- ncol(X)                                  # number of estimated coefficients
df_res <- n - p                               # residual degrees of freedom

#Coefficients: Estimate

XtX_inv <- solve(t(X) %*% X)
beta    <- XtX_inv %*% t(X) %*% y      # beta-hat = (X'X)^-1 X'y

# Residuals section
fitted_vals <- X %*% beta
res         <- y - fitted_vals
res_summary <- quantile(res, probs = c(0, .25, .5, .75, 1))   # Min 1Q Median 3Q Max

# Residual standard error
RSS    <- sum(res^2)
sigma2 <- RSS / df_res
sigma  <- sqrt(sigma2)                 # "Residual standard error"

# Coefficients: Std. Error, t value, Pr(>|t|)
se   <- sqrt(diag(sigma2 * XtX_inv))
tval <- beta / se
pval <- 2 * pt(abs(tval), df = df_res, lower.tail = FALSE)

coef_table <- cbind(Estimate = beta, `Std. Error` = se,
                    `t value` = tval, `Pr(>|t|)` = pval)

# R-squared, adjusted R-squared
has_int <- attr(terms(linear_model), "intercept") == 1
TSS <- if (has_int) sum((y - mean(y))^2) else sum(y^2)   # no-intercept case uses raw SS
R2  <- 1 - RSS / TSS
int <- as.integer(has_int)
adjR2 <- 1 - (1 - R2) * (n - int) / df_res

# F-statistic
df_model <- p - int
Fstat <- ((TSS - RSS) / df_model) / (RSS / df_res)
Fpval <- pf(Fstat, df_model, df_res, lower.tail = FALSE)

# Compare with summary() 
s <- summary(linear_model)

res_summary; summary(residuals(linear_model))          # quantiles (summary() also shows mean, ignore it)
coef_table;  s$coefficients
c(sigma, s$sigma)
c(R2, s$r.squared)
c(adjR2, s$adj.r.squared)
c(Fstat, s$fstatistic["value"])
c(Fpval, pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE))

# Diagnostics ----

# Perform visual diagnostics
par(mfrow = c(2, 2))   # Plot in a 2x2 grid
plot(LM, which = 1)    # Residuals vs fitted
qqPlot(LM, reps = 1e4) # QQ-plot
mtext("QQ-plot", 3, 0.25)
plot(LM, which = 3)    # Scale-location
plot(LM, which = 5)    # Cook's distance
par(mfrow = c(1, 1))   # Restore the default

## Residuals vs. fitted ----

plot(linear_model, which = 1)

linear_model |> 
  augment() |> 
  ggplot(aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey")

## Quantile-quantile plot ----

plot(linear_model, which = 2)

# Extraire les valeurs
valeurs <- df$height

# Ordonner les valeurs
valeurs <- sort(valeurs)

# Calcul des ranks
r <- rank(valeurs)

# Calcul des probabilités cumulées empiriques
p <- (r - 0.5) / length(valeurs)
p

# Calcul des quantiles normaux correspondant aux probabilités empiriques : 
# valeurs qui seraient observées si les données suivaient une loi normale
quantiles_normaux <- round(qnorm(p, mean = mean(valeurs), sd = sd(valeurs)), 1)
quantiles_normaux

# Calcul des quantiles normaux centrés réduits correspondant aux probabilités empiriques
quantiles_normaux_Std <- round(qnorm(p, mean = 0, sd = 1), 1)
quantiles_normaux_Std

# Assemblage dans un data frame
mydata <- data.frame(valeurs, r, p, quantiles_normaux, quantiles_normaux_Std)


ggplot(mydata, aes(x=quantiles_normaux_Std, y=valeurs))+
  geom_point()

shapiro.test(x = linear_model$residuals)

## Scale-location ----

plot(linear_model, which = 3)

linear_model |> 
  augment() |> 
  ggplot(aes(x = .fitted, y = sqrt(.std.resid))) +
  geom_point()

# Breusch-Pagan test to assess non-constant variance


lmtest::bptest(linear_model)

## Residuals vs. leverage ----

plot(linear_model, which = 5)

# Diagnostics using easystats ----

performance::check_model(linear_model)
performance::check_normality(linear_model)
performance::check_heteroscedasticity(linear_model)
