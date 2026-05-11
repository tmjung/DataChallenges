library(DBI)
library(RSQLite)

con <- dbConnect(SQLite(), "data/ffm_vfpa_eisenzeit.sqlite")

sql_text <- paste(readLines("data/ffm_vfpa_eisenzeit.sql", encoding = "UTF-8"), collapse = "\n")
stmts <- strsplit(sql_text, ";", fixed = TRUE)[[1]]
stmts <- trimws(stmts)
stmts <- stmts[nzchar(stmts) & !startsWith(stmts, "--")]

for (stmt in stmts) {
  dbExecute(con, stmt)
}

dbListTables(con)
dbDisconnect(con)


# Load the data into R
data <- read.csv("data/ffm_vfpa_eisenzeit.csv")

# Quick inspection
print(str(data))           # View structure

print("done")
head(data)          # See first rows
nrow(data)          # 15,395 sites
colnames(data)      # All 37 columns