find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)

  repeat {
    if (
      dir.exists(file.path(current, "src")) &&
        dir.exists(file.path(current, "data"))
    ) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Projektroot nicht gefunden. Bitte aus dem Projektordner oder src starten.")
    }

    current <- parent
  }
}

project_root <- find_project_root()
project_path <- function(...) file.path(project_root, ...)
