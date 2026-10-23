## =====================================================================
##  THE SIMULATION STUDY.  This is the only file you need to call.
##
##      source("run_simulation.R")
##      run_simulation("pilot")     # ~2 min: does it all wire up
##      run_simulation("bench")     # ~2 min: measure this machine's throughput
##      run_simulation("small")     # ~4 h:   400 subjects, all four scenarios
##      run_simulation("full")      # ~4 h:   the real 1,101-subject design
##
##      summarise_sim("small")      # re-draw tables and figures from saved replicates
##      sim_progress("small")       # from another R session, while it runs
##
##  Every replicate is checkpointed; a killed run resumes where it stopped.
##  See simulation/README.md for the protocol and the compute budget.
## =====================================================================
## Every path in the simulation is relative to the project folder, and the
## parallel workers re-source this file from there, so run it from that folder.
if (!file.exists("run_model.R") || !dir.exists("simulation"))
  stop("run this from the amyloid_mixture_cure project folder: setwd(\"/path/to/amyloid_mixture_cure\") first")
suppressMessages(source("run_model.R"))
for (f in sort(list.files("simulation", pattern = "\\.R$", full.names = TRUE))) source(f)
