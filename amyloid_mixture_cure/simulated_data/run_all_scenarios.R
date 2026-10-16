## =====================================================================
##  RUN EVERY SCENARIO EXCEPT BASELINE, one after another.
##
##  In RStudio, with the amyloid_mixture_cure project open, either
##    * paste this in the Console:
##        rstudioapi::jobRunScript("simulated_data/run_all_scenarios.R", workingDir = getwd())
##    * or open this file and choose  Source > Source as Background Job
##      (working directory: the project folder).
##  The Background Jobs tab then shows progress, and RStudio stays usable.
##
##  Each scenario fits 5 replicates in parallel (same as the baseline run),
##  about 1 hour per scenario. A replicate whose fit is already in
##  simulated_data/results is skipped, so if anything stops you can just
##  run this again and it picks up where it left off.
## =====================================================================

SCENARIOS <- c("sample_size_half", "onset_spread_lo", "onset_spread_hi",
               "speed_spread_lo", "speed_spread_hi", "meas_error_lo", "meas_error_hi")
N_SIM <- 5L
TIER  <- "infer"

source("simulated_data/fit_simulated_data.R")

## keep the Mac awake for as long as this R process runs (macOS caffeinate)
if (Sys.info()[["sysname"]] == "Darwin")
  system(sprintf("caffeinate -i -w %d", Sys.getpid()), wait = FALSE)

t_all <- Sys.time()
status <- data.frame(scenario = SCENARIOS, ok = NA, hours = NA_real_)
for (k in seq_along(SCENARIOS)) {
  s <- SCENARIOS[k]; t0 <- Sys.time()
  message(sprintf("\n==== %s (%d of %d), started %s ====", s, k, length(SCENARIOS),
                  format(t0, "%H:%M")))
  res <- try(run_simulation_study(n_sim = N_SIM, scenario = s, tier = TIER))
  status$ok[k] <- !inherits(res, "try-error")
  if (status$ok[k]) try(plot_sim_fits(scenario = s, tier = TIER))   # needs the qpdf package
  status$hours[k] <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
  message(sprintf("==== %s %s in %.1f h ====", s, if (status$ok[k]) "finished" else "FAILED",
                  status$hours[k]))
}

## the NMSE / NMISE tables for every scenario that now has fits (baseline included)
source("simulated_data/evaluate_simulation.R")
try(evaluate_all(tier = TIER))

message(sprintf("\nAll done in %.1f h", as.numeric(difftime(Sys.time(), t_all, units = "hours"))))
print(status, row.names = FALSE)
