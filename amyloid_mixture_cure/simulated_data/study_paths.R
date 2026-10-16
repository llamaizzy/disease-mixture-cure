## =====================================================================
##  WHERE THE SIMULATION STUDY KEEPS ITS FILES, and how a file records
##  what produced it.
##
##    simulated_data/data/<data_id>/            the simulated datasets
##        data_info.rds                         truth, scenarios, seed scheme, gen_hash
##        sim_<scenario>_repNN.rds
##        figures/                              sim_check_*.pdf
##    simulated_data/runs/<run_id>/             one model setup fitted to one data_id
##        run_info.rds                          data_id, created, git commit, notes
##        fits/<tier>_<scenario>_repNN_fit.rds  each carries fit$stamp
##        logs/   results/   figures/
##
##  A data_id changes only when the DATA-GENERATING setup changes (truth,
##  generator code, seed scheme). A run_id changes whenever the MODEL or its
##  settings change. Nothing is overwritten: a change is a new id.
##
##  STAMPS. Code is hashed after parsing, so comments and whitespace do not
##  count, and only the code that produces the file is included:
##    gen_hash   : R/00-02 + the generator functions in simulate_data.R + the truth
##    model_hash : R/00-04 + the fitting functions in fit_simulated_data.R + SETTINGS
##    data_hash  : the dataset's own contents
##  A fit is reused only if its stamp matches the current code; otherwise
##  run_simulation_study() stops rather than mixing code versions in one run.
## =====================================================================

SIM_ROOT <- "simulated_data"

data_dir <- function(data_id, ...) file.path(SIM_ROOT, "data", data_id, ...)
run_dir  <- function(run_id, ...)  file.path(SIM_ROOT, "runs", run_id, ...)

## ---- hashing ------------------------------------------------------------
.hash <- function(x) digest::digest(x, algo = "md5")

## parsed (comment-free) code of a file; with `names`, only those top-level assignments
.code_of <- function(file, names = NULL) {
  ex <- as.list(parse(file, keep.source = FALSE))
  if (!is.null(names)) {
    lhs <- vapply(ex, function(e) if (is.call(e) && as.character(e[[1]]) %in% c("<-", "="))
      paste(deparse(e[[2]]), collapse = "") else "", "")
    miss <- setdiff(names, lhs)
    if (length(miss)) stop(sprintf("%s does not define: %s", file, paste(miss, collapse = ", ")))
    ex <- ex[lhs %in% names]
  }
  vapply(ex, function(e) paste(deparse(e), collapse = "\n"), "")
}

.r_files <- function(upto) {
  f <- sort(list.files("R", "^0[0-9]_.*[.]R$", full.names = TRUE))
  f[as.integer(substr(basename(f), 1, 2)) <= upto]
}

GEN_FUNCTIONS <- c("get_baseline_truth", "make_scenarios", "check_truth_support", "rtnorm",
                   "get_design", "true_var_knots", "true_var_basis", "truth_curves",
                   "accum_levels", "draw_latent", "get_y", "get_dat", "true_x_tilde",
                   "simulate_one")
FIT_FUNCTIONS <- c("offgrid_at_inits", "subset_subjects", "configure_samplers",
                   "TIER_CFG", "fit_simulated_model")

gen_hash <- function(truth) .hash(list(
  code  = c(unlist(lapply(.r_files(2), .code_of)),
            .code_of(file.path(SIM_ROOT, "simulate_data.R"), GEN_FUNCTIONS)),
  truth = truth))

model_hash <- function(S = SETTINGS) .hash(list(
  code = c(unlist(lapply(.r_files(4), .code_of)),
           .code_of(file.path(SIM_ROOT, "fit_simulated_data.R"), FIT_FUNCTIONS)),
  settings = S))

data_hash <- function(sim) .hash(sim[c("dat", "truth", "latent", "seed")])

## git commit and whether the study code has uncommitted changes (recorded, not compared)
git_state <- function() {
  g <- function(...) tryCatch(suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
                              error = function(e) character(0))
  commit <- g("rev-parse", "HEAD")
  dirty  <- g("status", "--porcelain", "--", "R", file.path(SIM_ROOT, "*.R"))
  list(commit = if (length(commit)) commit[1] else NA_character_,
       dirty = length(dirty) > 0)
}

## ---- run bookkeeping ---------------------------------------------------
##  A run is tied to ONE data_id. The first call creates run_info.rds; later
##  calls must agree with it. Returns the run's data_id.
init_run <- function(run_id, data_id = NULL, notes = "") {
  if (missing(run_id) || !nzchar(run_id)) stop("name the run: run_id = \"<date>_<what changed>\"")
  f <- run_dir(run_id, "run_info.rds")
  if (file.exists(f)) {
    info <- readRDS(f)
    if (!is.null(data_id) && !identical(data_id, info$data_id))
      stop(sprintf("run '%s' is fitted to data '%s', not '%s'; use a new run_id",
                   run_id, info$data_id, data_id))
    return(info$data_id)
  }
  if (is.null(data_id)) stop(sprintf("run '%s' is new: say which data it uses (data_id = ...)", run_id))
  if (!file.exists(data_dir(data_id, "data_info.rds")))
    stop(sprintf("no data '%s': run simulate_all(data_id = \"%s\") first", data_id, data_id))
  for (d in c("fits", "logs", "results", "figures")) dir.create(run_dir(run_id, d), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(run_id = run_id, data_id = data_id, created = Sys.time(), git = git_state(),
               notes = notes), f)
  data_id
}

run_data_id <- function(run_id) {
  f <- run_dir(run_id, "run_info.rds")
  if (!file.exists(f)) stop(sprintf("no run '%s' (looked for %s)", run_id, f))
  readRDS(f)$data_id
}

## every run, with its data and notes
list_runs <- function() {
  ids <- basename(list.dirs(file.path(SIM_ROOT, "runs"), recursive = FALSE))
  do.call(rbind, lapply(ids, function(id) {
    i <- readRDS(run_dir(id, "run_info.rds"))
    data.frame(run_id = id, data_id = i$data_id, created = format(i$created, "%Y-%m-%d"),
               n_fits = length(list.files(run_dir(id, "fits"), "_fit[.]rds$")),
               notes = i$notes, stringsAsFactors = FALSE)
  }))
}

## ---- fit stamps ----------------------------------------------------------
STAMP_KEYS <- c("model_hash", "data_hash", "tier", "n_chains", "n_iter", "seed")

fit_stamp <- function(sim, tier, n_chains, n_iter, seed) {
  c(list(model_hash = model_hash(), data_hash = data_hash(sim), tier = tier,
         n_chains = as.integer(n_chains), n_iter = as.integer(n_iter), seed = as.integer(seed)),
    git = list(git_state()), R = R.version.string,
    nimble = as.character(utils::packageVersion("nimble")), time = Sys.time())
}

## which keys of an existing fit's stamp differ from `want` (character(0) = reusable)
stamp_mismatch <- function(stamp, want) {
  if (is.null(stamp)) return("no stamp (fit made before stamping)")
  STAMP_KEYS[!vapply(STAMP_KEYS, function(k) identical(stamp[[k]], want[[k]]), TRUE)]
}
