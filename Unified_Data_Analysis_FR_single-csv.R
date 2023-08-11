if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

max_temperature = 0
min_temperature = 0
lines_color <- c("green", "red", "blue", "orange", "purple")
i = 2

path <- "D:/Docs/RStudio/Test-R/data/spedizione/chimney/Oxic Pressure/Data10.csv"
data <- read.csv(path, row.names = NULL)
data$Survival = data$Alive / (data$Alive + data$Dead)
data_ordered <- data[order(data$Temperature), ]
xdata<-data[,1]
ydata<-data[,4]
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
  main = "LT50 Survival Curve"
)
points(xdata, ydata, col = lines_color[i])

start_p1 = 100
start_p2 = 30
start_p3 = 40
fit = nls(ydata ~ p1 / (1 + (xdata / p2) ^ p3), 
  start = list(p1 = start_p1, p2 = start_p2, p3 = start_p3),
  control = nls.control(maxiter = 500))
new = data.frame(xdata = seq(min(xdata), max(xdata), len = 200))
fit_params = summary(fit)
print(summary(fit))

lines(new$xdata, predict(fit, newdata = new), col = lines_color[i])