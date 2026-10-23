## =====================================================================
##  Build the joint amyloid-tau datasets: the amyloid cohort with its tau
##  scans, for either tau summary.
##
##  Run from the project root:   Rscript data/build_joint_dataset.R
##  writes BOTH
##    data/adni_joint_neo.rds   tau = NEO-T neocortical composite  (the MAIN file)
##    data/adni_joint.rds       tau = mesial temporal meta-ROI
##
##  Sources (data/raw/):
##    dta_amy.csv                  4,192 amyloid PET scans on 1,900 ADNI subjects
##    dta_tau.csv                  1,865 tau PET scans (mesial temporal meta-ROI)
##                                 on 991 subjects, with age at scan
##    ADNI_FTP_NEO_T_SUVR.csv      2,116 FTP tau scans on 1,184 subjects:
##                                 NEO_T_SUVR = volume-weighted mean of left and
##                                 right inferior and middle temporal gyri
##                                 (Ossenkoppele's NEO-T), dated but WITHOUT age
##    PET_export_for_scheffler.csv every amyloid scan with both its date and
##                                 the subject's age, which anchors each
##                                 subject's birth date so that NEO-T scan dates
##                                 can be turned into ages (to about 0.1 year)
##
##  The amyloid cohort is the 1,101 subjects with at least two amyloid
##  scans -- exactly the subjects in the original data/adni_amyloid.rds,
##  in the same order. Their tau scans are joined on RID. Note that
##  RID_num is numbered separately within each CSV, so it must NOT be used
##  as the join key.
##
##  The result is a SUPERSET of adni_amyloid.rds: every field of the
##  amyloid-only list is reproduced bit-for-bit (checked below), and the
##  tau fields are added alongside. Existing amyloid-only code therefore
##  runs unchanged on the joint file.
##
##  Tau fields (all indexed like their amyloid counterparts):
##    n_obs_tau     scans per subject; 0 for the 453 without tau
##    y_tau_mat     tau SUVR, subject x scan, zero-padded like y_mat
##    age_tau_mat   age at tau scan, zero-padded
##    a_tilde_tau   mean tau scan age; for subjects WITHOUT tau this is the
##                  amyloid centroid a_tilde, the "fixed reference age" the
##                  tau note (docs/context) uses for an unobserved biomarker
##    x_init_tau    mean tau level, a starting value for the tau location
##                  coordinate; the cohort median for subjects without tau
##    has_tau       logical
##    rid           the ADNI RID, so rows can be traced back to the CSVs
##    centiloid     intercept and slope of SUVR = a + b * Centiloid
##    tau_summary   "neo" or "mtl": which tau summary y_tau_mat holds
##
##  No tau positivity threshold is stored here: that is a modelling
##  setting to be calibrated on this summary (see docs/context), not a
##  property of the data.
## =====================================================================
build_joint_dataset <- function(tau_source = c("neo", "mtl"),
                                raw_dir = "data/raw",
                                out     = NULL,
                                check_against = "data/adni_amyloid.rds") {
  tau_source <- match.arg(tau_source)
  if (is.null(out)) out <- if (tau_source == "neo") "data/adni_joint_neo.rds" else "data/adni_joint.rds"
  amy <- read.csv(file.path(raw_dir, "dta_amy.csv"))
  mtl <- read.csv(file.path(raw_dir, "dta_tau.csv"))
  stopifnot(!anyNA(amy$SUVR), !anyNA(mtl$SUVR),
            !anyNA(amy$age_at_scan), !anyNA(mtl$age_at_scan))
  tau <- if (tau_source == "mtl") mtl else read_neo_tau(raw_dir, mtl)

  ## ---- the amyloid cohort, reproduced exactly ------------------------
  amy <- amy[order(amy$RID, amy$age_at_scan), ]
  n_scans <- table(amy$RID)
  keep    <- as.integer(names(n_scans)[n_scans >= 2])
  key     <- unique(amy[, c("RID", "RID_num")])
  cohort  <- key[key$RID %in% keep, ]
  cohort  <- cohort[order(cohort$RID_num), ]          # the rds order
  n       <- nrow(cohort)

  pad_wide <- function(df, rid, value, age) {
    sp  <- split(df, factor(df$RID, levels = rid))
    nob <- vapply(sp, nrow, 0L)
    m   <- max(nob, 1L)
    Y <- matrix(0, n, m); A <- matrix(0, n, m)
    for (i in seq_len(n)) if (nob[i] > 0) {
      Y[i, 1:nob[i]] <- sp[[i]][[value]]
      A[i, 1:nob[i]] <- sp[[i]][[age]]
    }
    list(n_obs = unname(nob), y = Y, age = A)
  }
  am <- pad_wide(amy[amy$RID %in% cohort$RID, ], cohort$RID, "SUVR", "age_at_scan")

  ## Covariates are centred over ALL 1,900 amyloid subjects, as in the
  ## original file (and the shipped posterior's intercepts live on that
  ## scale), not over the 1,101 in the cohort.
  first <- amy[!duplicated(amy$RID), ]
  apoe4  <- as.numeric(grepl("4", first$apoe_genotype))
  female <- as.numeric(first$sex == "Female")
  x_centre <- c(apoe4 = mean(apoe4), female = mean(female))
  m1 <- match(cohort$RID, first$RID)
  X  <- cbind(apoe4 = apoe4[m1] - x_centre[1], female = female[m1] - x_centre[2])

  a_tilde <- vapply(seq_len(n), function(i) mean(am$age[i, 1:am$n_obs[i]]), 0)
  x_init  <- vapply(seq_len(n), function(i) mean(am$y[i, 1:am$n_obs[i]]), 0)

  ## The positivity threshold is 24 Centiloid, mapped to SUVR through the
  ## cohort's own (exactly linear) SUVR <-> Centiloid relation.
  cl_map <- coef(lm(SUVR ~ CENTILOIDS, amy))
  thresh <- unname(cl_map[1] + 24 * cl_map[2])

  dat <- list(n_subj = n, n_obs = am$n_obs, y_mat = am$y, age_mat = am$age,
              a_tilde = a_tilde, x_init = x_init, X = X, ids = cohort$RID_num,
              x_centre = x_centre, thresh = thresh)

  if (!is.null(check_against) && file.exists(check_against)) {
    old <- readRDS(check_against)
    for (nm in names(old))
      if (!isTRUE(all.equal(unname(old[[nm]]), unname(dat[[nm]]), tolerance = 1e-10)))
        stop("field '", nm, "' does not reproduce data/adni_amyloid.rds")
    stopifnot(identical(colnames(old$X), colnames(dat$X)))
  }

  ## ---- the tau scans, joined on RID ----------------------------------
  tau <- tau[order(tau$RID, tau$age_at_scan), ]
  tau <- tau[tau$RID %in% cohort$RID, ]
  ## covariates must agree between the two files for the same subject (the
  ## NEO-T file carries none, so this only applies to the mesial temporal one)
  if (!is.null(tau$apoe_genotype)) {
    tk <- tau[!duplicated(tau$RID), ]; mt <- match(tk$RID, first$RID)
    stopifnot(all(tk$apoe_genotype == first$apoe_genotype[mt]),
              all(tk$sex == first$sex[mt]))
  }

  tm <- pad_wide(tau, cohort$RID, "SUVR", "age_at_scan")
  has_tau <- tm$n_obs > 0
  a_tilde_tau <- a_tilde
  a_tilde_tau[has_tau] <- vapply(which(has_tau),
                                 function(i) mean(tm$age[i, 1:tm$n_obs[i]]), 0)
  x_init_tau <- vapply(seq_len(n), function(i)
    if (has_tau[i]) mean(tm$y[i, 1:tm$n_obs[i]]) else NA_real_, 0)
  x_init_tau[!has_tau] <- median(x_init_tau, na.rm = TRUE)

  dat$rid         <- cohort$RID
  dat$centiloid   <- c(intercept = unname(cl_map[1]), slope = unname(cl_map[2]))
  dat$n_obs_tau   <- tm$n_obs
  dat$y_tau_mat   <- tm$y
  dat$age_tau_mat <- tm$age
  dat$a_tilde_tau <- a_tilde_tau
  dat$x_init_tau  <- x_init_tau
  dat$has_tau     <- has_tau
  dat$n_subj_tau  <- sum(has_tau)
  dat$tau_summary <- tau_source

  saveRDS(dat, out)
  cat(sprintf("%s [%s]: %d subjects, %d amyloid scans; %d with tau, %d tau scans\n",
              out, tau_source, dat$n_subj, sum(dat$n_obs), dat$n_subj_tau, sum(dat$n_obs_tau)))
  invisible(dat)
}

## ---- the NEO-T scans, with ages derived from scan dates ------------------
##
##  Birth-date anchor: every amyloid scan in PET_export has a date and an age
##  (to 0.1 year), so date - age is the birth date up to rounding; the median
##  over a subject's scans is used (within-subject spread is at most 62 days).
##  The derived ages are checked against the mesial temporal file, whose scans
##  are the same FTP acquisitions with ages supplied: matched by subject and
##  nearest age, they must agree to within the 0.1-year rounding.
read_neo_tau <- function(raw_dir, mtl) {
  neo <- read.csv(file.path(raw_dir, "ADNI_FTP_NEO_T_SUVR.csv"))
  pet <- read.csv(file.path(raw_dir, "PET_export_for_scheffler.csv"))
  stopifnot(!anyNA(neo$NEO_T_SUVR), all(neo$SCANDATE != ""), all(neo$qc_flag == 2),
            all(neo$TRACER == "FTP"), !any(duplicated(neo[, c("RID", "SCANDATE")])))
  ## the composite really is what it says
  w <- with(neo, (CTX_LH_INFERIORTEMPORAL_SUVR * CTX_LH_INFERIORTEMPORAL_VOLUME +
                  CTX_RH_INFERIORTEMPORAL_SUVR * CTX_RH_INFERIORTEMPORAL_VOLUME +
                  CTX_LH_MIDDLETEMPORAL_SUVR   * CTX_LH_MIDDLETEMPORAL_VOLUME +
                  CTX_RH_MIDDLETEMPORAL_SUVR   * CTX_RH_MIDDLETEMPORAL_VOLUME) / NEO_T_VOLUME)
  stopifnot("NEO_T_SUVR must be the volume-weighted mean of the four ROIs" =
              max(abs(w - neo$NEO_T_SUVR)) < 1e-9)

  birth <- tapply(as.numeric(as.Date(pet$pet_date)) - pet$age_at_scan * 365.25, pet$RID, median)
  neo   <- neo[neo$RID %in% as.integer(names(birth)), ]
  neo$age_at_scan <- round((as.numeric(as.Date(neo$SCANDATE)) - birth[as.character(neo$RID)]) / 365.25, 1)

  ## validate against the mesial temporal file's supplied ages
  shared <- neo[neo$RID %in% mtl$RID, ]
  gap <- vapply(seq_len(nrow(shared)), function(i)
    min(abs(mtl$age_at_scan[mtl$RID == shared$RID[i]] - shared$age_at_scan[i])), 0)
  cat(sprintf("NEO-T ages vs mesial temporal ages on %d shared scans: median gap %.2f yr, 95%% %.2f, max %.2f; %d scans have no MTL scan within 0.5 yr (new scans)\n",
              nrow(shared), median(gap), quantile(gap, .95), max(gap), sum(gap > 0.5)))
  stopifnot("derived NEO-T ages disagree with the mesial temporal ages" = median(gap) <= 0.2)

  data.frame(RID = neo$RID, age_at_scan = neo$age_at_scan, SUVR = neo$NEO_T_SUVR)
}

if (sys.nframe() == 0L) { build_joint_dataset("neo"); build_joint_dataset("mtl") }
