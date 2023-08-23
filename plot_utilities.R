if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

chart_subtitle_func <- function(dir_path) {
  chart_subtitle_parts <- unlist(strsplit(dir_path, "/"))
  chart_subtitle <- paste(
    chart_subtitle_parts[(length(chart_subtitle_parts) - 1):length(chart_subtitle_parts)],
    collapse = "/"
  )
  chart_subtitle <- gsub("/", " - ", chart_subtitle)
  chart_subtitle <- gsub("_", " ", chart_subtitle)
  chart_subtitle <- tools::toTitleCase(chart_subtitle)
  return(chart_subtitle)
}

save_plot_func <- function(plot, path, filename, width, height) {
  dpi <- 300
  width_in_inches <- width / dpi
  height_in_inches <- height / dpi
  # Current script directory
  script_path <- normalizePath(sys.frame(1)$ofile)
  script_dir <- dirname(script_path)
  plot_dir <- file.path(script_dir, path)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  ggsave(
    plot = plot,
    width = width_in_inches,
    height = height_in_inches,
    dpi = dpi,
    path = plot_dir,
    filename = filename,
    type = "cairo-png"
  )
}