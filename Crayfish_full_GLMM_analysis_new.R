# FULL PAIRED CRAYFISH ANALYSIS
# Open in RStudio and click Source. Select your final Excel workbook when asked.
# Install dependencies ONCE before running:
# install.packages(c("readxl", "dplyr", "tidyr", "glmmTMB", "DHARMa", "emmeans"))
# Outputs are saved in Crayfish_GLMM_results beside the Excel file.
# Original Excel data are never overwritten. No significance-based model selection.
# Script is supplied for local execution; it has not been run in this environment.

required <- c("readxl", "dplyr", "tidyr", "glmmTMB", "DHARMa", "emmeans")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(glmmTMB); library(DHARMa); library(emmeans)})

# SETTINGS -------------------------------------------------------------------
excel_file <- "" # Blank opens file chooser. Or use "C:/Users/Lenovo/Documents/Das et al.,RVS(5).xlsx".
output_directory <- "" # Blank creates output folder beside the workbook.
n_simulations <- 1000
random_SD_warning_threshold <- 0.01 # Practical flag; not a formal significance test.
set.seed(20261012)

run_analysis <- function() {
  if (!nzchar(excel_file)) excel_file <- file.choose()
  if (!file.exists(excel_file)) stop("Excel file not found.")
  if (!nzchar(output_directory)) output_directory <- file.path(dirname(excel_file), "Crayfish_GLMM_results")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  if (!"GLM_dataset" %in% excel_sheets(excel_file)) stop("Workbook does not contain GLM_dataset.")

  # IMPORT: first worksheet row is blank; next row holds headings.
  dat <- read_excel(excel_file, sheet = "GLM_dataset", skip = 1, .name_repair = "unique")
  names(dat) <- trimws(names(dat))
  needed <- c("ID", "Site", "Individual", "Organ", "Tissue_g", "MP_count", "Presence",
              "Fibers", "Fragments", "Films", "Foam", "Red", "Green", "Black", "Blue",
              "Yellow", "White", "Grey", "Small", "Medium", "Large")
  if (length(setdiff(needed, names(dat)))) stop("Missing columns: ", paste(setdiff(needed, names(dat)), collapse = ", "))
  dat <- dat %>% filter(!is.na(ID)) %>% select(all_of(needed))
  numeric_names <- setdiff(needed, c("Site", "Individual", "Organ"))
  dat <- dat %>% mutate(across(all_of(numeric_names), as.numeric),
                        across(c(Site, Individual, Organ), ~ trimws(as.character(.x))))
  dat$Organ[dat$Organ %in% c("Gills", "GILLS", "gills")] <- "Gill"
  if (anyNA(dat)) stop("Required observations contain missing/non-numeric values. Check records; do not replace them automatically with zeros.")
  if (any(dat$Tissue_g <= 0)) stop("Organ wet mass must be positive and in grams.")
  for (column in setdiff(numeric_names, "Tissue_g")) {
    if (any(dat[[column]] < 0 | dat[[column]] %% 1 != 0)) stop("Invalid integer data in ", column)
  }
  if (any(!dat$Site %in% c("S1", "S2", "S3")) || any(!dat$Organ %in% c("Gill", "GIT"))) stop("Unexpected site or organ label.")
  if (any(dat$Presence != as.integer(dat$MP_count > 0))) stop("Presence does not agree with particle count.")
  if (anyDuplicated(dat[c("Individual", "Organ")])) stop("Duplicate animal-organ records.")
  if (any(table(dat$Individual) != 2L)) stop("Each animal must have exactly two organ records.")
  if (any(vapply(split(dat$Site, dat$Individual), function(x) length(unique(x)), integer(1)) != 1L)) stop("An individual appears at more than one site.")
  prefix_site <- c(L = "S1", K = "S2", H = "S3")
  inferred <- unname(prefix_site[substr(dat$Individual, 1, 1)])
  if (anyNA(inferred) || any(dat$Site != inferred)) stop("Expected mapping is L=S1, K=S2, H=S3; check IDs and sites.")
  for (columns in list(c("Fibers", "Fragments", "Films", "Foam"),
                       c("Red", "Green", "Black", "Blue", "Yellow", "White", "Grey"),
                       c("Small", "Medium", "Large"))) {
    if (any(rowSums(dat[columns]) != dat$MP_count)) stop("Morphotype/color/size totals disagree with MP_count.")
  }
  dat <- dat %>% mutate(Site = factor(Site, levels = c("S1", "S2", "S3")),
                        Organ = factor(Organ, levels = c("Gill", "GIT")),
                        Individual = factor(Individual),
                        log_mass = log(Tissue_g), MP_per_g = MP_count/Tissue_g)
  write.csv(dat, file.path(output_directory, "Validated_organ_data.csv"), row.names = FALSE)
  write.csv(as.data.frame.matrix(table(dat$Site, dat$Organ)), file.path(output_directory, "Pairing_by_site.csv"))

  descriptive <- dat %>% group_by(Site, Organ) %>%
    summarise(n = n(), particle_total = sum(MP_count), positive = sum(Presence),
              occurrence_percent = 100*mean(Presence), mean_particles_g = mean(MP_per_g),
              sample_SD = sd(MP_per_g), median_particles_g = median(MP_per_g),
              pooled_particles_g = sum(MP_count)/sum(Tissue_g), .groups = "drop")
  individual_occurrence <- dat %>% group_by(Site, Individual) %>%
    summarise(positive = any(MP_count > 0), .groups = "drop") %>% group_by(Site) %>%
    summarise(animals = n(), positive = sum(positive), .groups = "drop") %>%
    mutate(percent = 100*positive/animals)
  write.csv(descriptive, file.path(output_directory, "Descriptive_statistics.csv"), row.names = FALSE)
  write.csv(individual_occurrence, file.path(output_directory, "Individual_occurrence.csv"), row.names = FALSE)

  sink(file.path(output_directory, "Model_results.txt"), split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("PAIRED CRAYFISH ANALYSIS\nInput:", normalizePath(excel_file, winslash = "/"), "\n")
  cat("Animals:", n_distinct(dat$Individual), " Organ samples:", nrow(dat), " Suspected particles:", sum(dat$MP_count), "\n")
  cat("Mass unit: grams. Counts are analysed as supplied, with no additional blank correction.\n")
  cat("Pairing check:\n"); print(table(table(dat$Individual))); print(table(dat$Site, dat$Organ))
  print(descriptive); print(individual_occurrence)

  # A failed fit is recorded; other models continue.
  ctrl <- glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
  fit_safe <- function(formula, family) {
    tryCatch(glmmTMB(formula, family = family, data = dat, control = ctrl),
             error = function(e) {cat("FIT FAILED:", conditionMessage(e), "\n"); NULL})
  }
  models <- list(
    NB2_mass = fit_safe(MP_count ~ Site * Organ + log_mass + (1 | Individual), nbinom2()),
    NB1_mass = fit_safe(MP_count ~ Site * Organ + log_mass + (1 | Individual), nbinom1()),
    NB2_offset = fit_safe(MP_count ~ Site * Organ + offset(log_mass) + (1 | Individual), nbinom2()),
    NB1_offset = fit_safe(MP_count ~ Site * Organ + offset(log_mass) + (1 | Individual), nbinom1()),
    Occurrence = fit_safe(Presence ~ Site * Organ + (1 | Individual), binomial()),
    NB2_mass_no_random = fit_safe(MP_count ~ Site * Organ + log_mass, nbinom2()),
    NB1_mass_no_random = fit_safe(MP_count ~ Site * Organ + log_mass, nbinom1()),
    Occurrence_no_random = fit_safe(Presence ~ Site * Organ, binomial())
  )
  # Explicit formulas avoid update() evaluation issues when this script is sourced
  # from inside run_analysis(). Each reduced model uses the same observations,
  # response, family, and random-effect structure as its corresponding full fit.
  reduced <- list(
    NB2_mass = list(
      additive = fit_safe(MP_count ~ Site + Organ + log_mass + (1 | Individual), nbinom2()),
      no_site = fit_safe(MP_count ~ Organ + log_mass + (1 | Individual), nbinom2())),
    NB1_mass = list(
      additive = fit_safe(MP_count ~ Site + Organ + log_mass + (1 | Individual), nbinom1()),
      no_site = fit_safe(MP_count ~ Organ + log_mass + (1 | Individual), nbinom1())),
    NB2_offset = list(
      additive = fit_safe(MP_count ~ Site + Organ + offset(log_mass) + (1 | Individual), nbinom2()),
      no_site = fit_safe(MP_count ~ Organ + offset(log_mass) + (1 | Individual), nbinom2())),
    NB1_offset = list(
      additive = fit_safe(MP_count ~ Site + Organ + offset(log_mass) + (1 | Individual), nbinom1()),
      no_site = fit_safe(MP_count ~ Organ + offset(log_mass) + (1 | Individual), nbinom1())),
    Occurrence = list(
      additive = fit_safe(Presence ~ Site + Organ + (1 | Individual), binomial()),
      no_site = fit_safe(Presence ~ Organ + (1 | Individual), binomial())),
    NB2_mass_no_random = list(
      additive = fit_safe(MP_count ~ Site + Organ + log_mass, nbinom2()),
      no_site = fit_safe(MP_count ~ Organ + log_mass, nbinom2())),
    NB1_mass_no_random = list(
      additive = fit_safe(MP_count ~ Site + Organ + log_mass, nbinom1()),
      no_site = fit_safe(MP_count ~ Organ + log_mass, nbinom1())),
    Occurrence_no_random = list(
      additive = fit_safe(Presence ~ Site + Organ, binomial()),
      no_site = fit_safe(Presence ~ Organ, binomial()))
  )
  valid_model <- function(m) {
    !is.null(m) && identical(as.integer(m$fit$convergence), 0L) && isTRUE(m$sdr$pdHess) && is.finite(as.numeric(logLik(m)))
  }
  get_random_SD <- function(m) {
    if (is.null(m)) return(NA_real_)
    vc <- VarCorr(m)$cond
    if (!length(vc)) return(NA_real_)
    value <- attr(vc[[1]], "stddev")
    if (is.null(value)) return(NA_real_)
    as.numeric(value[1])
  }
  status <- bind_rows(lapply(names(models), function(name) {
    m <- models[[name]]
    if (is.null(m)) return(data.frame(Model = name, usable = FALSE, convergence = NA_integer_,
      positive_definite_Hessian = FALSE, random_SD = NA_real_, near_zero_random_effect = NA,
      logLik = NA_real_, AIC = NA_real_))
    sd <- get_random_SD(m)
    data.frame(Model = name, usable = valid_model(m), convergence = m$fit$convergence,
      positive_definite_Hessian = isTRUE(m$sdr$pdHess), random_SD = sd,
      near_zero_random_effect = ifelse(is.na(sd), NA, sd < random_SD_warning_threshold),
      logLik = as.numeric(logLik(m)), AIC = AIC(m))
  }))
  cat("\nMODEL STATUS AND AIC\n"); print(status)
  cat("Compare AIC only among models of the same response and same observations.\n")
  cat("A near-zero random SD is a boundary warning, not a test of independence.\n")
  write.csv(status, file.path(output_directory, "Model_status_AIC.csv"), row.names = FALSE)
  reduced_status <- bind_rows(lapply(names(reduced), function(name) {
    bind_rows(lapply(names(reduced[[name]]), function(term) {
      m <- reduced[[name]][[term]]
      data.frame(Model = name, Reduced = term, usable = valid_model(m),
                 convergence = if (is.null(m)) NA_integer_ else m$fit$convergence,
                 positive_definite_Hessian = if (is.null(m)) FALSE else isTRUE(m$sdr$pdHess),
                 logLik = if (is.null(m)) NA_real_ else as.numeric(logLik(m)),
                 AIC = if (is.null(m)) NA_real_ else AIC(m))
    }))
  }))
  cat("\nREDUCED MODEL STATUS\n"); print(reduced_status)
  write.csv(reduced_status, file.path(output_directory, "Reduced_model_status.csv"), row.names = FALSE)

  lr_rows <- list()
  lr_test <- function(label, smaller, larger) {
    if (!valid_model(smaller) || !valid_model(larger)) {
      cat("Skipped LR test:", label, "because one model is unstable.\n"); return(invisible(NULL))
    }
    ll0 <- logLik(smaller); ll1 <- logLik(larger)
    df <- attr(ll1, "df") - attr(ll0, "df")
    statistic <- 2*(as.numeric(ll1)-as.numeric(ll0))
    if (df <= 0 || statistic < -1e-5) {cat("Invalid nested comparison:", label, "\n"); return(invisible(NULL))}
    statistic <- max(0, statistic)
    row <- data.frame(Test = label, LR_chisq = statistic, df = df, p_value = pchisq(statistic, df, lower.tail = FALSE))
    lr_rows[[length(lr_rows)+1L]] <<- row
    print(row)
  }

  for (name in names(models)) {
    m <- models[[name]]
    cat("\n\nMODEL:", name, "\n")
    if (is.null(m)) next
    print(summary(m)); print(VarCorr(m)); cat("Optimizer message:", m$fit$message, "\n")
    if (!valid_model(m)) {cat("UNSTABLE FIT: do not report its estimates or tests as final.\n"); next}
    co <- as.data.frame(summary(m)$coefficients$cond)
    co$term <- rownames(co); rownames(co) <- NULL
    co$CI_lower <- co$Estimate - qnorm(.975)*co$`Std. Error`
    co$CI_upper <- co$Estimate + qnorm(.975)*co$`Std. Error`
    write.csv(co, file.path(output_directory, paste0(name, "_coefficients.csv")), row.names = FALSE)

    # Full versus additive model tests the interaction jointly (two df).
    additive <- reduced[[name]]$additive
    lr_test(paste(name, "site-by-organ interaction"), additive, m)
    # Joint site differences across both organs, retaining mass and random structure.
    no_site <- reduced[[name]]$no_site
    lr_test(paste(name, "joint site effects including interaction"), no_site, m)

    # Conditional comparisons; no integration over random-effect distribution.
    at_values <- if ("log_mass" %in% all.vars(formula(m))) list(log_mass = mean(dat$log_mass)) else list()
    tryCatch({
      organ_emm <- emmeans(m, ~ Organ | Site, at = at_values)
      organ_contrasts <- as.data.frame(summary(pairs(organ_emm), infer = c(TRUE, TRUE),
                                               adjust = "none", type = "response"))
      # Holm across all three organ-within-site p-values; CIs remain unadjusted 95%.
      organ_contrasts$p_Holm <- p.adjust(organ_contrasts$p.value, method = "holm")
      cat("\nORGAN COMPARISONS: unadjusted 95% CIs; Holm-adjusted p-values across three sites\n")
      print(organ_contrasts)
      write.csv(organ_contrasts, file.path(output_directory, paste0(name, "_organ_comparisons.csv")), row.names = FALSE)
      site_emm <- emmeans(m, ~ Site | Organ, at = at_values)
      site_contrasts <- as.data.frame(summary(pairs(site_emm), infer = c(TRUE, TRUE),
                                              adjust = "none", type = "response"))
      site_contrasts$p_Holm <- p.adjust(site_contrasts$p.value, method = "holm")
      write.csv(site_contrasts, file.path(output_directory, paste0(name, "_site_comparisons.csv")), row.names = FALSE)
    }, error = function(e) cat("Comparison extraction failed:", conditionMessage(e), "\n"))
  }

  cat("\nMASS PROPORTIONALITY: offset fixes log-mass coefficient at one\n")
  lr_test("NB2: log-mass coefficient equals one", models$NB2_offset, models$NB2_mass)
  lr_test("NB1: log-mass coefficient equals one", models$NB1_offset, models$NB1_mass)
  if (length(lr_rows)) write.csv(bind_rows(lr_rows), file.path(output_directory, "Likelihood_ratio_tests.csv"), row.names = FALSE)
  cat("No ordinary chi-square LR test is used for random-effect variance: zero variance is a boundary.\n")

  # EXACT PAIRED OCCURRENCE CHECKS --------------------------------------------
  paired <- dat %>% select(Site, Individual, Organ, Presence) %>% pivot_wider(names_from = Organ, values_from = Presence)
  paired_tests <- bind_rows(lapply(levels(dat$Site), function(site) {
    z <- paired[paired$Site == site, ]
    gill_only <- sum(z$Gill == 1 & z$GIT == 0)
    git_only <- sum(z$Gill == 0 & z$GIT == 1)
    p <- if (gill_only + git_only > 0) binom.test(gill_only, gill_only + git_only, p = .5)$p.value else 1
    data.frame(Site = site, Gill_only = gill_only, GIT_only = git_only, p_exact = p)
  })) %>% mutate(p_Holm = p.adjust(p_exact, method = "holm"))
  cat("\nEXACT PAIRED OCCURRENCE TESTS\n"); print(paired_tests)
  write.csv(paired_tests, file.path(output_directory, "Paired_occurrence_tests.csv"), row.names = FALSE)

  # RESIDUAL DIAGNOSTICS ------------------------------------------------------
  diagnostic_rows <- list()
  for (name in c("NB2_mass", "NB1_mass", "NB2_offset", "NB1_offset", "Occurrence")) {
    m <- models[[name]]
    if (!valid_model(m)) next
    tryCatch({
      cat("\nSimulating residuals:", name, "\n")
      res <- simulateResiduals(fittedModel = m, n = n_simulations, refit = FALSE, seed = 20261012)
      saveRDS(res, file.path(output_directory, paste0(name, "_DHARMa.rds")))
      # Ensure graphics device is closed even if a plotting command fails.
      pdf(file.path(output_directory, paste0(name, "_residuals.pdf")), width = 8, height = 6)
      tryCatch({
        plot(res)
        plotResiduals(res, form = dat$log_mass)
        plotResiduals(res, form = dat$Site)
        plotResiduals(res, form = dat$Organ)
      }, finally = dev.off())
      tests <- list(Uniformity = testUniformity(res), Dispersion = testDispersion(res),
                    Zero_inflation = testZeroInflation(res))
      capture.output(lapply(tests, print), file = file.path(output_directory, paste0(name, "_diagnostics.txt")))
      for (test_name in names(tests)) {
        diagnostic_rows[[length(diagnostic_rows)+1L]] <- data.frame(Model = name,
          Test = test_name, p_value = tests[[test_name]]$p.value)
      }
    }, error = function(e) cat("Diagnostic error for", name, ":", conditionMessage(e), "\n"))
  }
  if (length(diagnostic_rows)) write.csv(bind_rows(diagnostic_rows), file.path(output_directory, "Diagnostic_tests.csv"), row.names = FALSE)
  saveRDS(list(full = models, reduced = reduced), file.path(output_directory, "Fitted_models.rds"))
  cat("\nSESSION INFORMATION\n"); print(sessionInfo())
  capture.output(sessionInfo(), file = file.path(output_directory, "R_session_info.txt"))
  writeLines(c(
    "READ THE MODEL STATUS TABLE FIRST. Convergence=0 and a positive-definite Hessian are basic numerical requirements, not proof of adequate fit.",
    "Count candidate: site*organ plus log organ mass and an animal random intercept. NB1/NB2 and offset fits are sensitivity comparisons. No model is selected automatically by p-value.",
    "Mass-covariate contrasts compare modelled counts at equal organ mass. Offset contrasts compare rates per gram under proportionality. Occurrence contrasts are odds ratios, not probability ratios.",
    "The organ-comparison CIs are unadjusted 95% intervals. The p_Holm column adjusts across all three organ comparisons; site comparisons are adjusted across their six comparisons.",
    "Joint site LR tests include site and site:organ terms; they are not identical to an averaged main-effect site test.",
    "Random-intercept variance near zero is a boundary warning. In particular, consider paired/clustered sensitivity analyses for occurrence rather than assuming observations become independent.",
    "A nonsignificant residual test does not prove model adequacy. Inspect plots and the biological assumptions. Do not add zero inflation solely because zeros occur.",
    "These are asymptotic likelihood-ratio and Wald tests. With 47 animals, uncertain/borderline effects may warrant a parametric-bootstrap check before submission.",
    "No normality test on raw particle counts is needed to justify a negative binomial GLMM. Report actual convergence, diagnostics, effect estimates and CIs.",
    "This script does not perform Raman, PHI/MPLI or particle-composition inference. Counts are visually suspected particles. Blank and recovery validation remain separate responsibilities.",
    "Send Model_results.txt, Model_status_AIC.csv, Reduced_model_status.csv, Likelihood_ratio_tests.csv, Diagnostic_tests.csv, and the residual PDFs for interpretation."
  ), file.path(output_directory, "READ_ME.txt"))
  cat("\nFINISHED. Output folder:", normalizePath(output_directory, winslash = "/"), "\n")
  invisible(list(data = dat, models = models, status = status, directory = output_directory))
}

result <- run_analysis()
dat <- result$data
models <- result$models
# Now you can inspect individual objects, for example:
# summary(models$NB2_mass)
# VarCorr(models$NB2_mass)
# result$status
