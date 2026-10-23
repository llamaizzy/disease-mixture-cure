## Shared loader for the 80k posterior analysis. Run scripts from the project root.
## Loads the model's R/ files, and ONLY the function definitions from
## run_model.R (it currently ends with interactive code that cannot be sourced).
for (f in sort(list.files("R", full.names = TRUE))) suppressMessages(source(f))
local({
  ex <- parse("run_model.R")
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                    is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
    eval(e, envir = globalenv())
})
suppressPackageStartupMessages({library(posterior); library(ggplot2); library(patchwork)})
OUT <- "posterior_80k"
