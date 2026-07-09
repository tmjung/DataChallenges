library(readr)

# -------------------------------------------------------------------------
# SQL dumps aus data/ in CSV-Dateien konvertieren
# -------------------------------------------------------------------------

input_files <- c(
  "data/bodenkarten_vfpa_eisenzeit.sql",
  "data/msqr_vfpa_eisenzeit.sql"
)

output_dir <- "data"
initialized_outputs <- new.env(parent = emptyenv())

parse_insert_header <- function(line) {
  match <- regexec(
    "^INSERT INTO `([^`]+)` \\((.*)\\) VALUES",
    line
  )
  parts <- regmatches(line, match)[[1]]

  if (length(parts) == 0) {
    return(NULL)
  }

  columns <- regmatches(parts[3], gregexpr("`[^`]+`", parts[3]))[[1]]
  columns <- gsub("`", "", columns, fixed = TRUE)

  list(
    table = parts[2],
    columns = columns
  )
}

split_sql_rows <- function(values_text) {
  rows <- character()
  row_start <- NA_integer_
  depth <- 0
  in_string <- FALSE
  escaped <- FALSE

  chars <- strsplit(values_text, "", fixed = TRUE)[[1]]

  for (i in seq_along(chars)) {
    ch <- chars[[i]]
    if (in_string) {
      if (escaped) {
        escaped <- FALSE
      } else if (ch == "\\") {
        escaped <- TRUE
      } else if (ch == "'") {
        in_string <- FALSE
      }
    } else {
      if (ch == "'") {
        in_string <- TRUE
      } else if (ch == "(") {
        if (depth == 0) {
          row_start <- i
        }
        depth <- depth + 1
      } else if (ch == ")") {
        depth <- depth - 1

        if (depth == 0) {
          rows <- c(rows, substring(values_text, row_start, i))
          row_start <- NA_integer_
        }
      }
    }
  }

  rows
}

split_sql_fields <- function(row_text) {
  row_text <- sub("^\\(", "", row_text)
  row_text <- sub("\\),?$", "", row_text)

  fields <- character()
  field_start <- 1L
  in_string <- FALSE
  escaped <- FALSE

  chars <- strsplit(row_text, "", fixed = TRUE)[[1]]

  for (i in seq_along(chars)) {
    ch <- chars[[i]]

    if (in_string) {
      if (escaped) {
        escaped <- FALSE
      } else if (ch == "\\") {
        escaped <- TRUE
      } else if (ch == "'") {
        in_string <- FALSE
      }
    } else {
      if (ch == "'") {
        in_string <- TRUE
      } else if (ch == ",") {
        fields <- c(fields, substring(row_text, field_start, i - 1L))
        field_start <- i + 1L
      }
    }
  }

  c(fields, substring(row_text, field_start))
}

clean_sql_value <- function(value) {
  value <- trimws(value)

  if (identical(toupper(value), "NULL")) {
    return(NA_character_)
  }

  if (grepl("^0x[0-9A-Fa-f]+$", value)) {
    return(value)
  }

  if (startsWith(value, "'") && endsWith(value, "'")) {
    value <- substring(value, 2, nchar(value) - 1)
    value <- gsub("\\\\'", "'", value)
    value <- gsub("\\\\n", "\n", value)
    value <- gsub("\\\\r", "\r", value)
    value <- gsub("\\\\t", "\t", value)
    value <- gsub("\\\\\\\\", "\\", value)
  }

  value
}

append_rows_to_csv <- function(rows, columns, output_file) {
  parsed_rows <- lapply(rows, function(row) {
    values <- split_sql_fields(row)
    values <- vapply(values, clean_sql_value, character(1))

    if (length(values) != length(columns)) {
      stop(
        sprintf(
          "Spaltenanzahl passt nicht: erwartet %d, gefunden %d in %s",
          length(columns),
          length(values),
          output_file
        )
      )
    }

    as.list(stats::setNames(values, columns))
  })

  df <- do.call(rbind.data.frame, c(parsed_rows, stringsAsFactors = FALSE))
  names(df) <- columns

  if (!exists(output_file, envir = initialized_outputs, inherits = FALSE)) {
    if (file.exists(output_file)) {
      file.remove(output_file)
    }

    assign(output_file, TRUE, envir = initialized_outputs)
  }

  write_csv(
    df,
    output_file,
    append = file.exists(output_file)
  )
}

convert_sql_file <- function(input_file, output_dir = "data") {
  cat("\nConverting:", input_file, "\n")

  con <- file(input_file, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  current_header <- NULL
  current_values <- character()
  created_files <- character()

  flush_insert <- function() {
    if (is.null(current_header) || length(current_values) == 0) {
      return(invisible(NULL))
    }

    values_text <- paste(current_values, collapse = "\n")
    values_text <- sub("^.* VALUES\\s*", "", values_text)
    values_text <- sub(";\\s*$", "", values_text)

    rows <- split_sql_rows(values_text)
    output_file <- file.path(output_dir, paste0(current_header$table, ".csv"))

    append_rows_to_csv(
      rows = rows,
      columns = current_header$columns,
      output_file = output_file
    )

    created_files <<- unique(c(created_files, output_file))
    invisible(NULL)
  }

  repeat {
    line <- readLines(con, n = 1, warn = FALSE)

    if (length(line) == 0) {
      flush_insert()
      break
    }

    header <- parse_insert_header(line)

    if (!is.null(header)) {
      flush_insert()
      current_header <- header
      current_values <- line
    } else if (!is.null(current_header)) {
      current_values <- c(current_values, line)

      if (grepl(";\\s*$", line)) {
        flush_insert()
        current_header <- NULL
        current_values <- character()
      }
    }
  }

  created_files
}

all_created_files <- unlist(
  lapply(input_files, convert_sql_file, output_dir = output_dir),
  use.names = FALSE
)

cat("\nDone.\n")
cat("Created CSV files:\n")
for (file in unique(all_created_files)) {
  cat(" -", file, "\n")
}
