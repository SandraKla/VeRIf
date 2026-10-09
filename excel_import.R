# Keep the normal readxl import behaviour and report incomplete headers clearly.
read_excel_dataset <- function(path, sheet) {
  dataset <- readxl::read_excel(
    path, sheet = sheet, .name_repair = "minimal"
  )
  labels <- trimws(names(dataset))
  if (length(labels) < 3L || anyNA(labels) || any(!nzchar(labels)) ||
      anyDuplicated(labels)) {
    stop(paste(
      "The first data row must contain at least three non-empty, unique column names.",
      "Merged cells may cause problems.",
      "If the column headers themselves are merged, unmerge them and give each column a name before uploading again."
    ), call. = FALSE)
  }
  names(dataset) <- labels
  as.data.frame(dataset, stringsAsFactors = FALSE)
}
