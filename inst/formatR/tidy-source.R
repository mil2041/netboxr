# use this script as part of a pre-commit hook to tidy up R source files
# e.g., add `Rscript -e 'source("inst/formatR/tidy-source.R")'`
#       to `.git/hooks/pre-commit`
library(rprojroot)
library(formatR)

#Sys.setlocale('LC_ALL','C')

options("formatR.arrow" = TRUE,
        "formatR.blank" = TRUE,
        "formatR.brace.newline" = FALSE,
        "formatR.comment" = TRUE,
        "formatR.indent" = 2)

files <- dir(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), pattern="[.][r|R]$", recursive=TRUE, full.names=TRUE)
lapply(files, function(f) {
  cat("FILE: ", f, "\n")
  result = tryCatch({
    formatR::tidy_source(f, arrow=TRUE, file=f, indent=2, width=80)
  }, error = function(e) {
    cat("ERROR: \n")
  })
})
