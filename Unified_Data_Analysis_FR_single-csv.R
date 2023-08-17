if (!require("stringr")) {
  install.packages("stringr")
  library(stringr)
}

if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

max_temperature <- 0
min_temperature <- 0
lines_color <- c("green", "red", "blue", "orange", "purple")
i <- 2

path <- "data/spedizione/chimney/Oxic Pressure/Data08.csv"
chart_subtitle_parts <- unlist(strsplit(path, "/"))
chart_subtitle <- paste(
  chart_subtitle_parts[(length(chart_subtitle_parts) - 2):length(chart_subtitle_parts)],
  collapse = "/"
)
#path <- "data/nioz/Data03.csv"
data <- read.csv(path, row.names = NULL)
data$Survival <- data$Alive / (data$Alive + data$Dead)
data_ordered <- data[order(data$Temperature), ]
xdata <- data[, 1]
ydata <- data[, 4]
if (min(data$Temperature) < min_temperature) {
  min_temperature <- min(data$Temperature)
}
if (max(data$Temperature) > max_temperature) {
  max_temperature <- max(data$Temperature)
}

# Create empty plot with correct axis values and labels
plot(1,
  type = "n",
  xlim = c(min_temperature, max_temperature),
  ylim = c(0.0, 1.0),
  xlab = "Temperature (°C)",
  ylab = "Proportional Survival",
  main = paste("LT50 Survival Curve\n", chart_subtitle)
)

start_p1 <- 0.5
start_p2 <- 33
start_p3 <- 77
fit <- nls(ydata ~ p1 / (1 + (xdata / p2)^p3),
  start = list(p1 = start_p1, p2 = start_p2, p3 = start_p3),
  control = nls.control(maxiter = 500)
)
new <- data.frame(xdata = seq(min(xdata), max(xdata), len = 200))
fit_params <- summary(fit)
cat("########## NON-LINEAR REGRESSION ##########\n\n")
print(summary(fit))

# Plot non-linear chart
model_dose <- fit_params$parameters[2]
rounded_dose <- round(model_dose, digits = 2)
model_dose_y <- predict(fit, newdata = data.frame(xdata = model_dose))
model_dose_stderr <- fit_params$parameters[2, 2]
lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])
points(jitter(xdata, factor = 0.5), ydata, col = lines_color[i])
text(
  x = rounded_dose,
  y = model_dose_y,
  labels = sprintf("%.2f°C", model_dose),
  pos = cos(pi * i) + 2,
  cex = 0.8
)
# Adding flex points
points(
  x = model_dose,
  y = model_dose_y,
  type = "p",
  pch = 19,
  col = lines_color[i]
)
# Adding std error ranges
arrows(
  x0 = model_dose - model_dose_stderr,
  y0 = model_dose_y,
  x1 = model_dose + model_dose_stderr,
  y1 = model_dose_y,
  code = 3,
  angle = 90,
  length = 0.05,
  col = lines_color[i]
)
