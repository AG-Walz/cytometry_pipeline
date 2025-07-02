library(ggplot2)
library(dplyr)
library(scales)

# Required Libraries: flowWorkspace, dplyr, openxlsx
# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("Functions.R")

x <- c(100, 1000, 10000, 100000, 340000,500000, 1272000)
y <- c(0.006, 0.03, 1, 60, 1800, 4500,39600)

# Prepare data for ggplot
df <- data.frame(x = x, y = y)

# LOESS fit for smooth curve
loess_fit <- loess(y ~ x, span = 0.75)
x_seq <- seq(min(x), max(x), length.out = 200)
y_loess <- predict(loess_fit, newdata = data.frame(x = x_seq))
df_loess <- data.frame(x = x_seq, y = y_loess)

ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 2) +
  geom_line(data = df_loess, aes(x = x, y = y), color = "red", size = 1) +
  scale_y_continuous(
    name = "Observed Computing Time",
    labels = time_formatter,
    breaks = c(0, 600, 1800, 3600, 18000, 36000)
  ) +
  scale_x_continuous(name = "Number of cells per data subset") +
  ggtitle("Fleiss' kappa computing time per cell number for 5 raters") +
  theme_minimal()

# Predict for a random x value
x_rand <- 1175000
y_pred_spline <- approx(x_seq, y_loess, xout = x_rand)$y
y_pred_loess <- predict(loess_fit, newdata = data.frame(x = x_rand))
cat(sprintf("Prediction for x = %d:\n  Spline: %.3f\n  LOESS: %.3f\n", x_rand, y_pred_spline, y_pred_loess))

summary(loess_fit)
