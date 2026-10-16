## =====================================================================
##  EVALUATE THE SIMULATION STUDY
##
##  Run from the package root, after fit_simulated_data.R has produced fits:
##
##      source("simulated_data/evaluate_simulation.R")
##      ev <- evaluate_scenario("baseline")          # one scenario
##      ev <- evaluate_all()                         # every scenario with fits
##
##
##  FIVE TABLES
##    Table 1  Subject-specific random effects: delta_i (acceleration factor),
##             alpha_i (departure age); centred NMSE, RMSE, Pearson, Spearman
##    Table 2  Variance components: sigma_alpha, sigma_delta, sigma_b
##    Table 3  Coefficients: phi (cure fraction), gamma (age of onset),
##             psi (rate of onset)
##    Table 4  Weibull nuisance: shape kappa = 1/sigma_alpha; scale
##             lambda_g = exp(x_g' gamma) for each of the four covariate groups
##    Table 5  Functional objects (NMISE): rate curve r_A(y);
##             measurement-error variance
##  METRICS:
##  Scalar parameter theta (Tables 2-4): NMSE = (1/R) sum_r error_r^2 / theta^2
##    CAUTION: dividing by theta^2 only means something for a positive SCALE
##    parameter (Tables 2 and 4). It explodes when |theta| is small and depends on an
##    arbitrary origin for an intercept. The coefficients (Table 3)
##    are therefore displayed as RMSE in their own units FIRST; their
##    NMSE follows, and the long file carries small_truth = TRUE on the rows
##    it cannot be read for.
##
##  CENTRED NMSE, subject-specific effect (Table 1), per replicate r:
##      NMSE_r = sum_i (est_i - true_i)^2 / sum_i (true_i - mean(true))^2,
##      NMSE = mean_r NMSE_r
##    normalised by the SPREAD of the true effects about their own mean, not
##    by their raw sum of squares: alpha_i sits around 60 years, so the raw
##    version is dominated by that location and changes with the origin age is
##    measured from. This means 0 is perfect and 1 is what giving every subject
##    the cohort mean would score (roughly 1 - R^2), so delta_i, alpha_i and
##    the scenarios are comparable. The denominator is still POOLED over
##    subjects: many true delta_i are near 0 and dividing subject by subject
##    would let a few dominate.
##    RMSE is in the effect's own units (delta_i: log speed; alpha_i: years).
##    Pearson / Spearman: correlation of est_i with true_i within a replicate,
##    averaged over replicates.
##    Only TRUE ACCUMULATORS (Z_i = 1) are scored.
##
##  NMISE, functional object f(y) (Table 5), per replicate r:
##      NISE_r = int_Y (f_hat_r(y) - f(y))^2 dy / int_Y f(y)^2 dy,  NMISE = mean_r NISE_r
##    f_hat_r is the pointwise posterior mean curve. Y is where the data inform
##    the curve (EVAL$r_range, EVAL$var_range). The measurement-error VARIANCE
##    at true level y is omega(y)^2 * df/(df - 2), Student-t with df = t_df.
##
##  coverage = share of replicates whose 95% credible interval holds the truth
##  mcse = Monte Carlo standard error of the average over replicates: sd_r(per-replicate value) / sqrt(R)
##        for RMSE, by the delta method. For a scalar's coverage this is
##        sqrt(p(1-p)/R) up to R/(R-1), so it reads 0 when every replicate covers.
##
##  OUTPUT
##    eval_<tier>_all_long.csv   ONE ROW PER SCENARIO x QUANTITY, every metric
##                               with its mcse. Everything else is built from it.
##    eval_<tier>_all_table<k>.csv   one table, one column per scenario, one
##                               block of rows per displayed metric (numeric).
##    The printed tables show "value (mcse)" per scenario;
##      *  the scenario changed this parameter's true value
##      +  worst R-hat over the replicates is above EVAL$rhat_flag
## =====================================================================

if (!exists("pool_paired")) for (f in sort(list.files("R", full.names = TRUE))) source(f)

EVAL <- list(
  estimator = "median",            # "mean" or "median" of the posterior
  burn_frac = 1/3,               # discarded from the start of EACH chain
  r_range   = c(NA, 1.30),       # SUVR range for the rate curve; NA = the true mu_b
  var_range = c(0.50, 1.30),     # SUVR range for the measurement-error variance
  small     = 0.15,              # a coefficient with |truth| below this is flagged small_truth
  rhat_flag = 1.1                # a printed cell is marked + when max R-hat exceeds this
)

TABLES <- c(table1 = "Table 1. Subject-specific random effects",
            table2 = "Table 2. Variance components",
            table3 = "Table 3. Coefficients",
            table4 = "Table 4. Weibull nuisance",
            table5 = "Table 5. Functional objects")

## Scenario label
SCEN <- c(baseline = "baseline", sample_size_half = "n_half",
          onset_spread_lo = "onset_lo", onset_spread_hi = "onset_hi",
          speed_spread_lo = "speed_lo", speed_spread_hi = "speed_hi",
          meas_error_lo = "meas_lo", meas_error_hi = "meas_hi")

## The metrics printed for each table (column of the long file = panel title).
##  bias is in the quantity's own units where the truth is fixed across
##  scenarios (Tables 1, 3); rel_bias = bias / |truth| where a scenario moves
##  the truth (Tables 2, 4), so the columns stay comparable.
PANELS <- list(
  table1 = c(NMSE = "NMSE, centred", rmse = "RMSE", bias = "bias",
             pearson = "Pearson correlation", spearman = "Spearman correlation",
             coverage = "95% interval coverage (share of subjects)"),
  table2 = c(NMSE = "NMSE", rel_bias = "relative bias, bias / truth",
             coverage = "95% interval coverage"),
  table3 = c(rmse = "RMSE", bias = "bias", coverage = "95% interval coverage",
             NMSE = "NMSE (unreliable where |truth| is small)"),
  table4 = c(NMSE = "NMSE", rel_bias = "relative bias, bias / truth",
             coverage = "95% interval coverage"),
  table5 = c(NMSE = "NMISE"))

## ---- small helpers ---------------------------------------------------
.est   <- function(x) if (EVAL$estimator == "median") median(x) else mean(x)
.trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
.burn  <- function(m, bf = EVAL$burn_frac) m[(floor(nrow(m) * bf) + 1):nrow(m), , drop = FALSE]
.mcse  <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
.ci    <- function(m) apply(m, 2, quantile, c(.025, .975))        # 95% interval, per column
.short <- function(s) ifelse(s %in% names(SCEN), SCEN[s], s)
.shown <- function(d) d[!grepl("departed before last scan", d$quantity), , drop = FALSE]  # kept in the files only
.sig   <- function(x) sub("\\.$", "", formatC(signif(x, 2), digits = 2, format = "fg", flag = "#"))

## R-hat on per-chain draws of one quantity (same formula as convergence())
.rhat <- function(L) {
  n <- min(lengths(L)); L <- lapply(L, function(x) x[seq_len(n)])
  M <- sapply(L, mean); V <- sapply(L, var)
  sqrt(((n - 1) / n * mean(V) + var(M)) / mean(V))
}

## ---- scalar quantities (Tables 2-4) -------------------------------------
##  Each is a function of the parameters, applied identically to the truth
##  and to every posterior draw, so a derived quantity (the Weibull shape and
##  scales) is estimated by the posterior mean of that quantity.
##
##  The Weibull departure-age model:
##      alpha_i = alpha_min + lambda_i * E_i^sigma_alpha,  E_i ~ Exp(1),
##      lambda_i = exp(x_i' gamma),  x_i = (1, APOE4 - c1, female - c2),
##  truncated to (alpha_min, alpha_max]. Its parameters beyond the gamma
##  coefficients are the SHAPE kappa = 1/sigma_alpha and the SCALE lambda_i.
##
##  The covariates are binary, so lambda_i takes four values, one per group.
##  alpha_min and alpha_max are fixed settings, not estimated.
GROUPS <- data.frame(label  = c("non-carrier male", "non-carrier female",
                                "APOE4 carrier male", "APOE4 carrier female"),
                     apoe4  = c(0, 0, 1, 1), female = c(0, 1, 0, 1))

scalar_spec <- function() {
  data.frame(
    item = c(rep("2", 3), rep("3", 8), rep("4", 5)),
    component = c(rep("variance", 3),
                  rep("cure fraction", 3), rep("age of onset", 3),
                  rep("rate of onset", 2),
                  "Weibull shape", rep("Weibull scale", 4)),
    quantity = c("sigma_alpha", "sigma_delta", "sigma_b",
                 "phi[1]", "phi[2]", "phi[3]",
                 "gamma[1]", "gamma[2]", "gamma[3]",
                 "psi[1]", "psi[2]",
                 "kappa = 1/sigma_alpha",
                 sprintf("lambda, %s", GROUPS$label)),
    stringsAsFactors = FALSE)
}

## m: matrix, one row per draw
scalar_values <- function(m, x_centre) {
  g <- function(nm) m[, nm]
  lam <- sapply(seq_len(nrow(GROUPS)), function(k)
    exp(g("gamma[1]") + g("gamma[2]") * (GROUPS$apoe4[k] - x_centre[1]) +
                        g("gamma[3]") * (GROUPS$female[k] - x_centre[2])))
  cbind(g("sigma_alpha"), g("sigma_delta"), g("sigma_b"),
        g("phi[1]"), g("phi[2]"), g("phi[3]"),
        g("gamma[1]"), g("gamma[2]"), g("gamma[3]"),
        g("psi[1]"), g("psi[2]"),
        1 / g("sigma_alpha"), matrix(lam, nrow = nrow(m)))
}

truth_matrix <- function(tr) {
  v <- c(setNames(tr$phi, sprintf("phi[%d]", 1:3)),
         setNames(tr$gamma, sprintf("gamma[%d]", 1:3)),
         setNames(tr$psi, sprintf("psi[%d]", 1:2)),
         sigma_alpha = tr$sigma_alpha, sigma_delta = tr$sigma_delta, sigma_b = tr$sigma_b)
  matrix(v, nrow = 1, dimnames = list(NULL, names(v)))
}

## ---- subject-level draws: alpha_i and delta_i (Table 1) -----------------
##  The fit monitors x_tilde_i and z_i, not alpha_i and delta_i; both are
##  rebuilt draw by draw exactly as subject_quantities() does (R/05_figures.R):
##      delta_i = x_i' psi + sigma_delta * z_i
##      alpha_i = a_tilde_i - exp(-delta_i) { G_ext(x_tilde_i) - G(mu_b) }
##  using the global draw ALIGNED with each thinned subject draw.
subject_draws <- function(fit) {
  bA <- fit$bA; S <- fit$settings; dat <- fit$dat; n <- dat$n_subj
  pp <- pool_paired(fit, burn_frac = EVAL$burn_frac); g <- pp$g; u <- pp$u
  XT <- u[, sprintf("x_tilde[%d]", 1:n), drop = FALSE]
  ZZ <- u[, sprintf("z[%d]", 1:n), drop = FALSE]
  K  <- nrow(g); AL <- DE <- matrix(NA_real_, K, n)
  for (k in seq_len(K)) {
    rgr <- rate_grid(bA, g[k, sprintf("theta[%d]", 1:bA$K)])
    G   <- clock_grid(bA, rgr, S$thresh)
    mb  <- g[k, "mu_b"]; Gb <- sa_interp(bA$y, G, mb)
    rb  <- max(sa_interp(bA$y, rgr, mb), S$r_min)
    delta  <- as.numeric(dat$X %*% g[k, grep("^psi", colnames(g))]) + g[k, "sigma_delta"] * ZZ[k, ]
    g0  <- ifelse(XT[k, ] >= mb, sa_interp(bA$y, G, XT[k, ]), Gb + (XT[k, ] - mb) / rb)
    DE[k, ] <- delta
    AL[k, ] <- dat$a_tilde - exp(-delta) * (g0 - Gb)
  }
  list(alpha = AL, delta = DE)
}

## ---- one fit against its truth -----------------------------------------
##  Returns the per-replicate rows: scalar (Tables 2-4), subject (Table 1),
##  functional (Table 5). Every row carries error / sq_error (where the
##  quantity has units), nse (the normalised squared error) and covered.
evaluate_one <- function(fit) {
  S <- fit$settings; tr <- fit$truth; dat <- fit$dat
  id <- data.frame(scenario = fit$scenario, rep = fit$rep, stringsAsFactors = FALSE)
  burned <- lapply(fit$chains, .burn)
  rhat_of <- function(nms) max(sapply(nms, function(nm) .rhat(lapply(burned, function(m) m[, nm]))))

  ## Tables 2-4: scalars --------------------------------------------------
  sp <- scalar_spec()
  xc <- as.numeric(dat$x_centre)
  per_chain <- lapply(burned, scalar_values, x_centre = xc)
  pooled <- do.call(rbind, per_chain)
  est <- apply(pooled, 2, .est); q <- .ci(pooled)
  tv  <- as.numeric(scalar_values(truth_matrix(tr), xc))
  rh  <- sapply(seq_along(tv), function(j) .rhat(lapply(per_chain, function(m) m[, j])))
  scal <- cbind(id, sp, truth = tv, estimate = est, error = est - tv, sq_error = (est - tv)^2,
                nse = (est - tv)^2 / tv^2, covered = as.numeric(tv >= q[1, ] & tv <= q[2, ]), rhat = rh)

  ## Table 1: subject-specific random effects -------------------------------
  sd_  <- subject_draws(fit)
  acc  <- which(fit$latent$Z == 1)
  last <- dat$age_mat[cbind(seq_len(dat$n_subj), dat$n_obs)]
  seen <- acc[fit$latent$alpha[acc] < last[acc]]
  subj <- list()
  for (v in c("delta", "alpha")) {
    e <- apply(sd_[[v]], 2, .est); q <- .ci(sd_[[v]]); tru <- fit$latent[[v]]
    lab <- if (v == "delta") "delta_i (acceleration factor)" else "alpha_i (departure age)"
    for (who in list(list(i = acc,  tag = "all true accumulators"),
                     list(i = seen, tag = "departed before last scan"))) {
      i <- who$i; err <- e[i] - tru[i]
      subj[[length(subj) + 1]] <- cbind(id, data.frame(
        item = "1", component = lab, quantity = sprintf("%s: %s", sub(" .*", "", lab), who$tag),
        n_subjects = length(i),
        error = mean(err), sq_error = mean(err^2),
        nse = sum(err^2) / sum((tru[i] - mean(tru[i]))^2),          # CENTRED
        pearson = cor(e[i], tru[i]), spearman = cor(e[i], tru[i], method = "spearman"),
        covered = mean(tru[i] >= q[1, i] & tru[i] <= q[2, i]),
        stringsAsFactors = FALSE))
    }
  }
  subj <- do.call(rbind, subj)

  ## Table 5: functional objects ----------------------------------------
  draws <- do.call(rbind, burned)
  bA <- fit$bA; vA <- fit$vA
  th <- sprintf("theta[%d]", 1:bA$K); nu <- sprintf("nu[%d]", 1:vA$K)
  lo <- if (is.na(EVAL$r_range[1])) tr$mu_b else EVAL$r_range[1]
  ky <- bA$y >= lo & bA$y <= EVAL$r_range[2]
  Bk <- bA$B[ky, , drop = FALSE]
  r_true <- as.numeric(exp(Bk %*% tr$theta))
  Rd     <- exp(Bk %*% t(draws[, th]))
  r_hat  <- if (EVAL$estimator == "median") apply(Rd, 1, median) else rowMeans(Rd)
  infl <- S$t_df / (S$t_df - 2)
  kv <- vA$y >= EVAL$var_range[1] & vA$y <= EVAL$var_range[2]
  Bv <- vA$B[kv, , drop = FALSE]
  v_true <- infl * as.numeric(exp(2 * Bv %*% tr$nu))
  Vd     <- infl * exp(2 * Bv %*% t(draws[, nu]))
  v_hat  <- if (EVAL$estimator == "median") apply(Vd, 1, median) else rowMeans(Vd)
  func <- cbind(id, data.frame(
    item = "5", component = "functional object",
    object = c("rate curve r_A(y)", "measurement-error variance omega(y)^2 df/(df-2)"),
    y_range = c(sprintf("[%.3f, %.2f]", lo, EVAL$r_range[2]),
                sprintf("[%.2f, %.2f]", EVAL$var_range[1], EVAL$var_range[2])),
    nise = c(.trapz(bA$y[ky], (r_hat - r_true)^2) / .trapz(bA$y[ky], r_true^2),
             .trapz(vA$y[kv], (v_hat - v_true)^2) / .trapz(vA$y[kv], v_true^2)),
    rhat = c(rhat_of(th[-1]), rhat_of(nu)),          # worst coefficient of the curve; theta[1] is fixed
    stringsAsFactors = FALSE))

  list(scalar = scal, subject = subj, functional = func)
}

## ---- summarise the replicates: one row per quantity ----------------------
##  $long   : every quantity of every table, identical columns (NA where a
##            metric does not apply). For table5 the NMSE column is the NMISE.
##  $tables : the same rows split by table, without the empty columns.
summarise_tables <- function(reps) {
  scal <- do.call(rbind, lapply(reps, `[[`, "scalar"))
  subj <- do.call(rbind, lapply(reps, `[[`, "subject"))
  func <- do.call(rbind, lapply(reps, `[[`, "functional"))
  func$quantity <- sprintf("%s on %s", sub("variance .*", "variance", func$object), func$y_range)
  func$nse <- func$nise

  one <- function(d) {
    col <- function(nm) if (nm %in% names(d)) d[[nm]] else rep(NA_real_, nrow(d))
    mse <- mean(col("sq_error"))
    data.frame(
      scenario = .short(d$scenario[1]), table = paste0("table", d$item[1]),
      component = d$component[1], quantity = d$quantity[1],
      n_sim = nrow(d), n_subjects = mean(col("n_subjects")),
      truth = col("truth")[1], mean_estimate = mean(col("estimate")),
      bias = mean(col("error")), bias_mcse = .mcse(col("error")),
      rel_bias = mean(col("error")) / abs(col("truth")[1]),
      rel_bias_mcse = .mcse(col("error")) / abs(col("truth")[1]),
      rmse = sqrt(mse), rmse_mcse = .mcse(col("sq_error")) / (2 * sqrt(mse)),
      NMSE = mean(d$nse), NMSE_mcse = .mcse(d$nse),
      pearson = mean(col("pearson")), pearson_mcse = .mcse(col("pearson")),
      spearman = mean(col("spearman")), spearman_mcse = .mcse(col("spearman")),
      coverage = mean(col("covered")), coverage_mcse = .mcse(col("covered")),
      max_rhat = max(col("rhat")),
      stringsAsFactors = FALSE)
  }
  long <- do.call(rbind, lapply(list(subj, scal, func), function(d)
    do.call(rbind, lapply(split(d, factor(d$quantity, unique(d$quantity))), one))))
  long$small_truth <- long$table == "table3" & abs(long$truth) < EVAL$small
  rownames(long) <- NULL

  tables <- lapply(split(long, long$table), function(x) {
    x <- x[, !vapply(x, function(v) all(is.na(v)), TRUE), drop = FALSE]
    if (x$table[1] == "table5") names(x) <- sub("^NMSE", "NMISE", names(x))
    if (x$table[1] != "table3") x$small_truth <- NULL
    x$table <- NULL; rownames(x) <- NULL; x
  })
  list(long = long, tables = tables,
       by_sim = list(scalar = scal, subject = subj, functional = func))
}

## one scenario: every column of each table
print_tables <- function(tabs, header = "") {
  op <- options(width = 200); on.exit(options(op))
  if (nzchar(header)) cat("\n", header, "\n", sep = "")
  for (k in names(tabs)) {
    cat("\n", TABLES[[k]], "\n", sep = "")
    x <- .shown(tabs[[k]]); x$scenario <- NULL; x$component <- NULL; x$n_sim <- NULL
    print(format(x, digits = 3), row.names = FALSE)
  }
  invisible(tabs)
}

## ---- one scenario --------------------------------------------------------
##  Writes eval_<tier>_<scenario>_table1.csv ... _table5.csv,
##  and the per-replicate values to eval_<tier>_<scenario>_by_sim_*.csv.
evaluate_scenario <- function(scenario = "baseline", tier = "infer",
                              out_dir = "simulated_data/results", write = TRUE, verbose = TRUE) {
  files <- list.files(out_dir, sprintf("^%s_%s_rep[0-9]+_fit\\.rds$", tier, scenario), full.names = TRUE)
  if (!length(files)) { if (verbose) message("no fits for ", tier, " / ", scenario); return(NULL) }
  reps <- lapply(files, function(f) {
    if (verbose) message("  scoring ", basename(f))
    evaluate_one(readRDS(f))
  })
  res <- summarise_tables(reps)
  if (write) {
    for (k in names(res$tables))
      write.csv(res$tables[[k]], file.path(out_dir, sprintf("eval_%s_%s_%s.csv", tier, scenario, k)),
                row.names = FALSE)
    for (nm in names(res$by_sim))
      write.csv(res$by_sim[[nm]], file.path(out_dir, sprintf("eval_%s_%s_by_sim_%s.csv", tier, scenario, nm)),
                row.names = FALSE)
  }
  if (verbose) print_tables(res$tables, sprintf("%s / %s: %d replicates (estimator: posterior %s)",
                                                scenario, tier, length(files), EVAL$estimator))
  invisible(res)
}

## ---- across scenarios ------------------------------------------------------
## rows = quantity, columns = scenario, for one column of the long file
.wide <- function(d, v, lev) {
  q <- unique(d$quantity)
  m <- matrix(NA_real_, length(q), length(lev), dimnames = list(q, lev))
  m[cbind(match(d$quantity, q), match(d$scenario, lev))] <- as.numeric(d[[v]])
  m
}

## The printed tables: one block per metric, "value (mcse)" per scenario,
## * where the scenario changed the truth, + where R-hat is above the flag.
print_across <- function(long) {
  op <- options(width = 250); on.exit(options(op))
  lev <- unique(long$scenario)
  cat(sprintf("\nCells are value (mcse).  * the scenario changed this true value   + max R-hat > %s\n",
              EVAL$rhat_flag))
  for (k in names(PANELS)) {
    d <- .shown(long[long$table == k, ])
    if (!nrow(d)) next
    star <- .wide(d, "truth_changed", lev); star <- !is.na(star) & star == 1
    dag  <- .wide(d, "max_rhat", lev);      dag  <- !is.na(dag) & dag > EVAL$rhat_flag
    for (v in names(PANELS[[k]])) {
      est <- .wide(d, v, lev); se <- .wide(d, paste0(v, "_mcse"), lev)
      cell <- ifelse(is.na(est), "",
                     paste0(.sig(est), ifelse(is.na(se), "", sprintf(" (%s)", .sig(se))),
                            ifelse(star, "*", ""), ifelse(dag, "+", "")))
      out <- data.frame(quantity = rownames(est), stringsAsFactors = FALSE)
      if (k == "table3") out$truth <- .sig(d$truth[match(out$quantity, d$quantity)])
      out <- cbind(out, matrix(cell, nrow(est), dimnames = list(NULL, lev)), stringsAsFactors = FALSE)
      cat(sprintf("\n%s: %s\n", TABLES[[k]], PANELS[[k]][[v]]))
      print(format(out, justify = "left"), row.names = FALSE, right = FALSE)
    }
  }
  invisible(long)
}

## ---- every scenario that has fits ------------------------------------------
##  Writes eval_<tier>_all_long.csv (one row per scenario x quantity) and, built
##  from it, eval_<tier>_all_table1.csv ... _table5.csv (one column per scenario,
##  one block of rows per displayed metric).
evaluate_all <- function(tier = "infer", scenarios = NULL, out_dir = "simulated_data/results") {
  if (is.null(scenarios)) {
    f <- list.files(out_dir, sprintf("^%s_.+_rep[0-9]+_fit\\.rds$", tier))
    scenarios <- unique(sub(sprintf("^%s_(.+)_rep[0-9]+_fit\\.rds$", tier), "\\1", f))
  }
  scenarios <- c(intersect(names(SCEN), scenarios), setdiff(scenarios, names(SCEN)))   # display order
  res <- setNames(lapply(scenarios, evaluate_scenario, tier = tier, out_dir = out_dir,
                         verbose = FALSE), scenarios)
  res <- res[!vapply(res, is.null, TRUE)]
  if (!length(res)) stop("no fits found in ", out_dir)

  long <- do.call(rbind, lapply(res, `[[`, "long")); rownames(long) <- NULL
  lev  <- unique(long$scenario)
  ## did the scenario move this parameter's true value away from the baseline's?
  key  <- paste(long$table, long$quantity); b <- long$scenario == "baseline"
  bt   <- long$truth[b][match(key, key[b])]
  long$truth_changed <- !is.na(bt) & !is.na(long$truth) & abs(long$truth - bt) > 1e-8 * abs(bt)
  write.csv(long, file.path(out_dir, sprintf("eval_%s_all_long.csv", tier)), row.names = FALSE)

  wide <- lapply(setNames(names(PANELS), names(PANELS)), function(k) {
    d <- .shown(long[long$table == k, ])
    w <- do.call(rbind, lapply(names(PANELS[[k]]), function(v) {
      m <- .wide(d, v, lev)
      data.frame(metric = if (k == "table5") "NMISE" else v, quantity = rownames(m), m,
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    rownames(w) <- NULL
    write.csv(w, file.path(out_dir, sprintf("eval_%s_all_%s.csv", tier, k)), row.names = FALSE)
    w
  })
  print_across(long)
  invisible(list(by_scenario = res, long = long, wide = wide))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  evaluate_all(tier = if (length(args) >= 1) args[1] else "infer")
}
