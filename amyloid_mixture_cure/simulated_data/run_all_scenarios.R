## =====================================================================
##  RUN EVERY SCENARIO, one after another, into ONE run.
##
##  In RStudio, with the amyloid_mixture_cure project open, either
##    * paste this in the Console:
##        rstudioapi::jobRunScript("simulated_data/run_all_scenarios.R", workingDir = getwd())
##    * or open this file and choose  Source > Source as Background Job
##      (working directory: the project folder).
##  The Background Jobs tab then shows progress, and RStudio stays usable.
##
##  Each scenario fits N_SIM replicates in parallel. Fits go to
##  simulated_data/runs/<RUN_ID>/; the datasets come from
##  simulated_data/data/<DATA_ID>/ and are made there if missing. A replicate
##  whose fit is already in the run (and matches the current code) is skipped,
##  so if anything stops you can just run this again and it picks up where it
##  left off. A CHANGE to the model or settings needs a new RUN_ID; a change to
##  how the data are simulated needs a new DATA_ID.
## =====================================================================

RUN_ID  <- "2026-10-16_crn-20reps"   # "<date>_<what changed>"; a new name for every change
DATA_ID <- "crn-v1"                  # common random numbers across scenarios (new); "v1" is the old data
NOTES   <- "common seeds across scenarios, 20 replicates, model as in 2026-10-07_v1"
SCENARIOS <- "baseline"
# SCENARIOS <- c("baseline", "sample_size_half", "onset_spread_lo", "onset_spread_hi",
#               "speed_spread_lo", "speed_spread_hi", "meas_error_lo", "meas_error_hi")
N_SIM <- 20L
TIER  <- "infer"

source("simulated_data/fit_simulated_data.R")

## keep the Mac awake for as long as this R process runs (macOS caffeinate)
if (Sys.info()[["sysname"]] == "Darwin")
  system(sprintf("caffeinate -i -w %d", Sys.getpid()), wait = FALSE)

## the datasets first (a no-op for any already on disk); a new DATA_ID uses common seeds
simulate_all(n_rep = N_SIM, which = SCENARIOS, data_id = DATA_ID, verbose = FALSE)

t_all <- Sys.time()
status <- data.frame(scenario = SCENARIOS, ok = NA, hours = NA_real_)
for (k in seq_along(SCENARIOS)) {
  s <- SCENARIOS[k]; t0 <- Sys.time()
  message(sprintf("\n==== %s (%d of %d), started %s ====", s, k, length(SCENARIOS),
                  format(t0, "%H:%M")))
  res <- try(run_simulation_study(n_sim = N_SIM, scenario = s, tier = TIER,
                                  run_id = RUN_ID, data_id = DATA_ID, notes = NOTES))
  status$ok[k] <- !inherits(res, "try-error")
  if (status$ok[k]) try(plot_sim_fits(scenario = s, tier = TIER, run_id = RUN_ID))   # needs the qpdf package
  status$hours[k] <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  message(sprintf("==== %s %s in %.1f h ====", s, if (status$ok[k]) "finished" else "FAILED",
                  status$hours[k]))
}

## the tables for every scenario that now has fits in this run
try(evaluate_all(tier = TIER, run_id = RUN_ID))

message(sprintf("\nAll done in %.1f h", as.numeric(difftime(Sys.time(), t_all, units = "hours"))))
print(status, row.names = FALSE)
