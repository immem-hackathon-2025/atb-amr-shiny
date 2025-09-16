options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))
if (!requireNamespace("duckdb", quietly = TRUE)) {
  install.packages("duckdb")
}
