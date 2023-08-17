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