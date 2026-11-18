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

# t-value = (estimate - 0) / standard error

(summary(linear_model)$coefficients[1] - 0) / summary(linear_model)$coefficients[3]
summary(linear_model)$coefficients[5]

(summary(linear_model)$coefficients[2] - 0) / summary(linear_model)$coefficients[4]
summary(linear_model)$coefficients[6]

pt(q = summary(linear_model)$coefficients[5], df = 29, lower.tail = F)


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
