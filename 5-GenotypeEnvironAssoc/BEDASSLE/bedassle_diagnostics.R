# Post-processing shared by the runner and the summary-only mode.
# Diagnostics are flags for human review, not a declaration of convergence.

validate_bedassle_settings <- function(ngen, burnin, samplefreq, savefreq,
                                      min_saved_states = 50000) {
  vals <- c(ngen, burnin, samplefreq, savefreq, min_saved_states)
  if (any(!is.finite(vals)) || any(vals != floor(vals)) ||
      ngen <= 0 || samplefreq <= 0 || savefreq <= 0 || burnin < 0 ||
      min_saved_states < 0 || burnin >= ngen) {
    stop("Invalid iteration, burn-in, sampling, checkpoint, or retained-state settings.")
  }
  if (ngen %% samplefreq != 0 || ngen %% savefreq != 0 ||
      savefreq %% samplefreq != 0) {
    stop("ngen must be divisible by samplefreq and savefreq; savefreq must be divisible by samplefreq.")
  }
  retained <- ngen / samplefreq - floor(burnin / samplefreq)
  if (retained < max(100, min_saved_states)) {
    stop(sprintf("Only %d post-burn-in states would remain; require at least %d.",
                 retained, max(100, min_saved_states)))
  }
  invisible(retained)
}

load_bedassle_chain <- function(path, burnin) {
  if (!file.exists(path)) stop("Missing MCMC output: ", path)
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  required <- c("ngen", "samplefreq", "a0", "aD", "aE", "a2", "beta", "phi_mat", "last.params")
  if (!all(required %in% ls(e))) stop("Incomplete BEDASSLE object: ", path)
  p <- e$last.params
  counter_names <- c("a0_moves", "aD_moves", "aE_moves", "a2_moves",
                     "beta_moves", "thetas_moves", "mu_moves", "phi_moves", "k", "loci")
  if (!all(counter_names %in% names(p))) stop("Cannot verify checkpoint completion: ", path)
  # Each iteration updates exactly one parameter family. Theta/mu counters
  # increment once per locus and phi counters once per population.
  updates <- p$a0_moves + p$aD_moves + p$aE_moves + p$a2_moves + p$beta_moves +
    p$thetas_moves / p$loci + p$mu_moves / p$loci + p$phi_moves / p$k
  if (!is.finite(updates) || abs(updates - (e$ngen - 1)) > 0.01) {
    stop("MCMC output is an unfinished checkpoint, not a completed chain: ", path)
  }
  if (e$ngen %% e$samplefreq != 0 || length(e$aD) != e$ngen / e$samplefreq) {
    stop("Saved-state dimensions do not match the recorded iterations: ", path)
  }
  n <- length(e$aD)
  aE <- if (is.null(dim(e$aE))) matrix(e$aE, nrow = 1) else e$aE
  if (nrow(aE) != 1 || ncol(aE) != n) stop("Expected one environmental predictor per model.")
  if (any(vapply(list(e$a0, e$a2, e$beta), length, integer(1)) != n) ||
      ncol(e$phi_mat) != n || nrow(e$phi_mat) != p$k) stop("Inconsistent saved parameter dimensions.")
  pops <- rownames(p$counts)
  if (is.null(pops)) pops <- seq_len(p$k)
  phi <- t(e$phi_mat)
  colnames(phi) <- paste0("phi_", pops)
  values <- cbind(a0 = e$a0, aD = e$aD, aE = as.numeric(aE),
                  a2 = e$a2, beta = e$beta, phi)
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("Saved chain contains non-finite or non-positive model parameters: ", path)
  }
  values <- cbind(values, aE_aD = values[, "aE"] / values[, "aD"])
  if (any(!is.finite(values))) stop("Non-finite environmental/geographic ratios: ", path)
  # BEDASSLE stores state i at iteration i * samplefreq.
  iterations <- seq_len(n) * e$samplefreq
  keep <- which(iterations > burnin)
  if (!is.finite(burnin) || burnin < 0 || burnin >= e$ngen || length(keep) < 100) {
    stop("Burn-in must leave at least 100 saved states: ", path)
  }
  chain <- coda::mcmc(values[keep, , drop = FALSE],
                      start = iterations[keep[1]], thin = e$samplefreq)
  list(chain = chain, full_values = values, iterations = iterations,
       ngen = e$ngen, samplefreq = e$samplefreq, burnin = burnin,
       retained = length(keep), last.params = p)
}

bedassle_parameter_diagnostics <- function(chain, ess_min = 400) {
  mat <- as.matrix(chain)
  rows <- lapply(colnames(mat), function(parameter) {
    x <- mat[, parameter]
    mc <- coda::mcmc(x)
    ess <- tryCatch(as.numeric(coda::effectiveSize(mc)), error = function(e) NA_real_)
    geweke <- tryCatch(as.numeric(coda::geweke.diag(mc)$z), error = function(e) NA_real_)
    heidel <- tryCatch(as.numeric(coda::heidel.diag(mc)[1, ]),
                       error = function(e) rep(NA_real_, 6))
    autocorrelation <- function(lag) {
      if (sd(x) == 0) return(NA_real_)
      as.numeric(stats::acf(x, plot = FALSE, lag.max = lag)$acf[lag + 1])
    }
    q <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
    flags <- character()
    if (!is.finite(ess) || ess < ess_min) flags <- c(flags, "low_or_unavailable_ESS")
    if (!is.finite(geweke) || abs(geweke) > 1.96) flags <- c(flags, "Geweke_flag")
    if (!is.finite(heidel[1]) || heidel[1] != 1) flags <- c(flags, "stationarity_flag")
    if (is.finite(heidel[2]) && heidel[2] > 1) flags <- c(flags, "additional_discard_suggested")
    if (!is.finite(heidel[4]) || heidel[4] != 1) flags <- c(flags, "halfwidth_flag")
    data.frame(Parameter = parameter, Mean = mean(x), Median = q[2],
               CrI_lower = q[1], CrI_upper = q[3],
               ESS = ess, MCSE_mean = if (is.finite(ess) && ess > 0) sd(x) / sqrt(ess) else NA_real_,
               ACF_lag1 = autocorrelation(1), ACF_lag10 = autocorrelation(10),
               ACF_lag50 = autocorrelation(50), Geweke_Z = geweke,
               Heidel_stationarity_pass = heidel[1], Heidel_start_draw = heidel[2],
               Heidel_p = heidel[3], Heidel_halfwidth_pass = heidel[4],
               Diagnostic_status = if (length(flags)) "review_required" else "no_numeric_flags",
               Diagnostic_flags = paste(flags, collapse = ";"), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

summarize_bedassle_file <- function(path, variable, burnin, diagnostic_prefix,
                                    ess_min = 400, min_saved_states = 50000,
                                    make_plots = TRUE) {
  obj <- load_bedassle_chain(path, burnin)
  dir.create(dirname(diagnostic_prefix), recursive = TRUE, showWarnings = FALSE)
  diagnostics <- bedassle_parameter_diagnostics(obj$chain, ess_min)
  diagnostics$Variable <- variable
  diagnostics$Retained_states <- obj$retained
  diagnostics$Burnin_iterations <- burnin
  write.csv(diagnostics, paste0(diagnostic_prefix, "_parameter_diagnostics.csv"), row.names = FALSE)
  write.csv(diagnostics[, c("Parameter", "ESS")], paste0(diagnostic_prefix, "_ess.csv"), row.names = FALSE)
  p <- obj$last.params
  families <- c("aD", "aE", "a2", "phi", "thetas", "mu")
  acceptance <- data.frame(Parameter = families,
    Proposals = vapply(families, function(x) p[[paste0(x, "_moves")]], numeric(1)),
    Accepted = vapply(families, function(x) p[[paste0(x, "_accept")]], numeric(1)))
  acceptance$Acceptance_rate <- with(acceptance, ifelse(Proposals > 0, Accepted / Proposals, NA_real_))
  write.csv(acceptance, paste0(diagnostic_prefix, "_acceptance.csv"), row.names = FALSE)
  if (make_plots) {
    # Explicitly disable ggs' automatic burn-in treatment: already excluded above.
    plot_data <- ggmcmc::ggs(obj$chain, burnin = FALSE)
    # ggs(burnin=FALSE) resets the burn-in label. Restore iteration offsets
    # for report axes only, without removing any additional saved states.
    attr(plot_data, 'nBurnin') <- start(obj$chain) - coda::thin(obj$chain)
    ggmcmc::ggmcmc(plot_data,
      file = paste0(diagnostic_prefix, "_diagnostics.pdf"), param_page = 5,
      plot = c("density", "traceplot", "running", "autocorrelation", "geweke"),
      simplify_traceplot = NULL)
    # Also retain a full-chain view so the specified burn-in boundary is visible.
    grDevices::pdf(paste0(diagnostic_prefix, "_full_chain_traces.pdf"), width = 10, height = 8)
    tryCatch({
      par(mfrow = c(3, 2), mar = c(3, 4, 2, 1))
      for (parameter in colnames(obj$full_values)) {
        plot(obj$iterations, obj$full_values[, parameter], type = "l",
             xlab = "MCMC iteration", ylab = parameter, main = variable)
        abline(v = burnin, col = "red", lty = 2)
      }
    }, finally = grDevices::dev.off())
  }
  ratio <- diagnostics[diagnostics$Parameter == "aE_aD", ]
  enough <- obj$retained >= min_saved_states
  status <- if (!enough || any(diagnostics$Diagnostic_status == "review_required"))
    "review_required" else "no_numeric_flags"
  summary <- data.frame(Variable = variable, aE_aD_Mean = ratio$Mean,
    aE_aD_Median = ratio$Median, CrI_lower = ratio$CrI_lower, CrI_upper = ratio$CrI_upper,
    CI = sprintf("[%.6g - %.6g]", ratio$CrI_lower, ratio$CrI_upper),
    CI_type = "95% equal-tailed posterior credible interval",
    Iterations = obj$ngen, Burnin_iterations = burnin, Samplefreq = obj$samplefreq,
    Retained_states = obj$retained, Retained_state_target_met = enough,
    Ratio_ESS = ratio$ESS, Min_parameter_ESS = min(diagnostics$ESS),
    Parameters_flagged = sum(diagnostics$Diagnostic_status == "review_required"),
    Diagnostic_status = status, Chains = 1L, Rhat = NA_real_, stringsAsFactors = FALSE)
  if (status == "review_required") warning("Review diagnostics before interpreting model: ", variable)
  list(summary = summary, diagnostics = diagnostics)
}
