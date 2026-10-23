## =====================================================================
##  The driver: tiers, the worker pool, checkpointing, and the budget.
##
##  THE BUDGET, measured on an 11-core M3 Pro (one process):
##      compile once            ~70 s (full cohort), ~35 s (400 subjects)
##      one iteration           1.06e-4 s per subject  (0.117 s at n = 1101)
##  Replicates run in parallel, one compiled model per worker, so the wall
##  clock is roughly  (replicates x chains x iterations x s/iter) / workers
##  times a contention factor -- about 1.3 with 8 workers, since the M3 Pro
##  has 5 performance and 6 efficiency cores and the efficiency cores are
##  slower. Run tier "bench" to measure that factor on your own machine; the
##  planner uses it if it exists.
##
##  Every replicate is saved as soon as it is scored. Rerunning a tier picks
##  up where it stopped, so a run you have to kill costs you at most one
##  replicate per worker.
## =====================================================================
suppressPackageStartupMessages(library(parallel))

SIM_SEC_PER_ITER_PER_SUBJ <- 1.06e-4          # measured; see header
SIM_CONTENTION_DEFAULT    <- 1.3

## n_rep is per scenario, in scenario order.
sim_tiers <- function() list(
  ## does everything wire up, end to end, in a few minutes
  pilot = list(n_sub = 250L, n_chains = 2L, n_iter = 400L,
               scenarios = c("fitted", "rate_only", "low_cure", "null"), n_rep = c(2L, 1L, 1L, 1L)),
  ## measures this machine's throughput with every worker busy
  bench = list(n_sub = 400L, n_chains = 2L, n_iter = 200L,
               scenarios = "fitted", n_rep = NA_integer_),
  ## the study at a reduced cohort: all four scenarios, an afternoon
  small = list(n_sub = 400L, n_chains = 3L, n_iter = 3000L,
               scenarios = c("fitted", "rate_only", "low_cure", "null"), n_rep = c(100L, 40L, 40L, 40L)),
  ## the study at the real design, the primary scenario only: overnight
  full  = list(n_sub = NA_integer_, n_chains = 3L, n_iter = 6000L,
               scenarios = "fitted", n_rep = 40L))

sim_plan <- function(cfg, n_workers, dir = "results/sim") {
  n_sub <- if (is.na(cfg$n_sub)) 1101L else cfg$n_sub
  bench <- file.path(dir, "bench_rate.rds")
  spi <- if (file.exists(bench)) readRDS(bench)$sec_per_iter_per_subj * 1.0 else
           SIM_SEC_PER_ITER_PER_SUBJ * SIM_CONTENTION_DEFAULT
  per_rep <- cfg$n_chains * cfg$n_iter * n_sub * spi / 60 + 0.4     # + scoring
  total <- sum(cfg$n_rep) * per_rep / n_workers + 1.5                  # + compile
  list(n_sub = n_sub, per_rep_min = per_rep, total_min = total, measured = file.exists(bench))
}

## ---- one job, run inside a worker that holds compiled model W ---------
sim_job <- function(job, cfg, W, dir) {
  f  <- file.path(dir, sprintf("rep_%s_%03d.rds", job$scenario, job$rep))
  fe <- file.path(dir, sprintf("err_%s_%03d.rds", job$scenario, job$rep))
  if (file.exists(f)) return("skip")
  t0 <- Sys.time()
  res <- tryCatch({
    truth <- sim_truth(job$scenario)
    sim   <- simulate_cohort(W$dat, truth, W$bA, W$vA, seed = job$seed)
    fit   <- sim_fit(W, sim$dat, cfg$n_chains, cfg$n_iter, seed = job$seed + 500000L)
    sc    <- score_replicate(fit, sim, truth)
    sc$scenario <- job$scenario; sc$rep <- job$rep; sc$seed <- job$seed
    sc$cfg <- cfg[c("n_sub", "n_chains", "n_iter")]
    saveRDS(sc, f); if (file.exists(fe)) unlink(fe)
    sprintf("done  %.1f min  worst R-hat %.2f", sc$minutes, max(sc$rhat, na.rm = TRUE))
  }, error = function(e) {
    saveRDS(list(error = TRUE, message = conditionMessage(e), job = job), fe)
    paste("ERROR", conditionMessage(e))
  })
  cat(sprintf("[%s] %-9s rep %3d  %s\n", format(Sys.time(), "%H:%M:%S"), job$scenario, job$rep, res))
  flush.console()
  res
}

## ---- the entry point --------------------------------------------------
run_simulation <- function(tier = c("pilot", "bench", "small", "full"),
                           scenarios = NULL, n_rep = NULL, n_iter = NULL, n_chains = NULL,
                           n_workers = NULL, dir = "results/sim", dry_run = FALSE,
                           summarise = TRUE) {
  tier <- match.arg(tier)
  cfg  <- sim_tiers()[[tier]]
  if (!is.null(scenarios)) { cfg$scenarios <- match.arg(scenarios, SIM_SCENARIOS, several.ok = TRUE)
                             if (is.null(n_rep)) cfg$n_rep <- rep(cfg$n_rep[1], length(scenarios)) }
  if (!is.null(n_rep))    cfg$n_rep    <- rep_len(as.integer(n_rep), length(cfg$scenarios))
  if (!is.null(n_iter))   cfg$n_iter   <- as.integer(n_iter)
  if (!is.null(n_chains)) cfg$n_chains <- as.integer(n_chains)
  if (is.null(n_workers)) n_workers <- max(1L, min(8L, detectCores() - 3L))
  if (tier == "bench") cfg$n_rep <- n_workers          # one job per worker, all at once
  out <- file.path(dir, tier); dir.create(out, showWarnings = FALSE, recursive = TRUE)

  jobs <- do.call(rbind, lapply(seq_along(cfg$scenarios), function(k)
    data.frame(scenario = cfg$scenarios[k], rep = seq_len(cfg$n_rep[k]),
               seed = 10000L * k + seq_len(cfg$n_rep[k]), stringsAsFactors = FALSE)))
  done <- file.exists(file.path(out, sprintf("rep_%s_%03d.rds", jobs$scenario, jobs$rep)))
  plan <- sim_plan(cfg, n_workers, dir)
  cat(sprintf("tier %s: %d subjects, %d chains x %d iterations, %d workers\n", tier, plan$n_sub,
              cfg$n_chains, cfg$n_iter, n_workers))
  for (k in seq_along(cfg$scenarios)) cat(sprintf("  %-10s %3d replicates\n", cfg$scenarios[k], cfg$n_rep[k]))
  cat(sprintf("  %d of %d already done; estimated %.0f min per replicate per worker, ~%.1f h wall clock for the rest (%s)\n",
              sum(done), nrow(jobs), plan$per_rep_min,
              sum(!done) * plan$per_rep_min / n_workers / 60 + 0.025, 
              if (plan$measured) "from your bench run" else "from the M3 Pro constant; run tier 'bench' to calibrate"))
  cat(sprintf("  progress: tail -f %s\n", file.path(out, "log.txt")))
  if (dry_run) return(invisible(jobs))
  jobs <- jobs[!done, , drop = FALSE]
  root <- normalizePath(".")

  if (nrow(jobs)) {
    t0 <- Sys.time()
    if (n_workers == 1L) {
      W <- sim_compile(sim_design(cfg$n_sub), make_rate_basis(), sim_var_basis())
      for (j in seq_len(nrow(jobs))) sim_job(jobs[j, ], cfg, W, out)
    } else {
      cl <- makeCluster(n_workers, outfile = file.path(out, "log.txt"))
      on.exit(stopCluster(cl), add = TRUE)
      clusterExport(cl, c("cfg", "root", "out"), envir = environment())
      cat("compiling on every worker ...\n"); flush.console()
      clusterEvalQ(cl, {
        setwd(root); suppressMessages(source("run_simulation.R"))
        nimbleOptions(verbose = FALSE, showCompilerOutput = FALSE, MCMCprogressBar = FALSE)
        W <- sim_compile(sim_design(cfg$n_sub), make_rate_basis(), sim_var_basis())
        cat(sprintf("[%s] worker %d compiled\n", format(Sys.time(), "%H:%M:%S"), Sys.getpid())); NULL })
      cat(sprintf("compiled (%.1f min); running %d replicates ...\n",
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(jobs))); flush.console()
      clusterApplyLB(cl, split(jobs, seq_len(nrow(jobs))),
                     function(job) { r <- sim_job(job, cfg, W, out); gc(); r })
    }
    cat(sprintf("finished in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  if (tier == "bench") {
    reps <- load_replicates(tier, dir); n_sub <- plan$n_sub
    spi <- vapply(reps, function(r) r$minutes * 60 / (r$n_iter * r$n_chains) / n_sub, 0)
    saveRDS(list(sec_per_iter_per_subj = mean(spi), n_workers = n_workers, date = Sys.time()),
            file.path(dir, "bench_rate.rds"))
    cat(sprintf("bench: %.2e s per iteration per subject with %d workers busy (%.2fx the single-process rate). Saved.\n",
                mean(spi), n_workers, mean(spi) / SIM_SEC_PER_ITER_PER_SUBJ))
    return(invisible(spi))
  }
  errs <- list.files(out, pattern = "^err_", full.names = TRUE)
  if (length(errs)) { cat(sprintf("%d replicates FAILED (rerun the tier to retry):\n", length(errs)))
    for (e in errs) cat("  ", basename(e), ": ", readRDS(e)$message, "\n") }
  if (summarise) summarise_sim(tier, dir) else invisible(NULL)
}

## How far along is a tier that is running in another process?
sim_progress <- function(tier, dir = "results/sim") {
  out <- file.path(dir, tier); fs <- list.files(out, pattern = "^rep_.*\\.rds$", full.names = TRUE)
  cat(sprintf("%d replicates finished", length(fs)))
  if (length(fs) > 1) {
    mt <- sort(file.mtime(fs)); rate <- as.numeric(difftime(max(mt), min(mt), units = "mins")) / (length(fs) - 1)
    cat(sprintf(", one every %.1f min recently", rate))
  }
  cat("\n"); invisible(length(fs))
}
