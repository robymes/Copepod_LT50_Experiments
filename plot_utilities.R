if (!require("MASS")) {
  install.packages("MASS")
  library(MASS)
}

if (!require("plotly")) {
  install.packages("plotly")
  library(plotly)
}

chart_subtitle_func <- function(dir_path) {
  chart_subtitle_parts <- unlist(strsplit(dir_path, "/"))
  chart_subtitle <- paste(
    chart_subtitle_parts[(length(chart_subtitle_parts) - 1):length(chart_subtitle_parts)],
    collapse = "/"
  )
  chart_subtitle <- gsub("/", " / ", chart_subtitle)
  chart_subtitle <- gsub("_", " ", chart_subtitle)
  chart_subtitle <- tools::toTitleCase(chart_subtitle)
  return(chart_subtitle)
}

save_plot_func <- function(plot, path, filename, width, height) {
  dpi <- 300
  width_in_inches <- width / dpi
  height_in_inches <- height / dpi
  ggsave("grafico.png",
    plot,
    width = width_in_inches,
    height = height_in_inches,
    dpi = dpi,
    path = path,
    filename = filename,
    antialias = "cairo"
  )
}