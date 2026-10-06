#!/usr/bin/env Rscript
# Deterministic tests for the scientific post-processing, without long MCMC runs.
arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
test_dir <- dirname(normalizePath(sub('^--file=', '', arg[1])))
source(file.path(dirname(test_dir), 'bedassle_diagnostics.R'))
stopifnot(requireNamespace('coda', quietly = TRUE))
temp <- tempfile('bedassle-tests-')
dir.create(temp)
expect_error <- function(expr) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'))
}
make_fixture <- function(path, incomplete = FALSE, corrupt = FALSE) {
  ngen <- 1000L; samplefreq <- 2L
  a0 <- aD <- a2 <- beta <- rep(1, 500)
  aE <- matrix(seq_len(500), nrow = 1)
  phi_mat <- matrix(1, nrow = 2, ncol = 500)
  if (corrupt) aD[500] <- 0
  last.params <- list(k = 2L, loci = 3L, counts = matrix(1, 2, 3),
    a0_moves = 0, aD_moves = if (incomplete) 998 else 999,
    aE_moves = 0, a2_moves = 0, beta_moves = 0, thetas_moves = 0,
    mu_moves = 0, phi_moves = 0, aD_accept = 100, aE_accept = 0,
    a2_accept = 0, thetas_accept = 0, mu_accept = 0, phi_accept = 0)
  save(ngen, samplefreq, a0, aD, aE, a2, beta, phi_mat, last.params, file = path)
}
tryCatch({
  stopifnot(validate_bedassle_settings(17000000, 4250000, 250, 100000) == 51000)
  expect_error(validate_bedassle_settings(2000000, 500000, 250, 100000))
  expect_error(validate_bedassle_settings(17000001, 4250000, 250, 100000))
  expect_error(validate_bedassle_settings(1000, 1000, 2, 100, 0))
  path <- file.path(temp, 'fixture.Robj')
  make_fixture(path)
  obj <- load_bedassle_chain(path, 200)
  stopifnot(obj$retained == 400, start(obj$chain) == 202,
            load_bedassle_chain(path, 201)$retained == 400)
  result <- suppressWarnings(summarize_bedassle_file(path, 'fixture', 200,
    file.path(temp, 'fixture'), min_saved_states = 100, make_plots = FALSE))
  s <- result$summary
  # Known quantiles of retained values 101:500. These are distribution
  # intervals, which must remain much wider than a confidence interval for a mean.
  stopifnot(abs(s$aE_aD_Mean - 300.5) < 1e-10,
            abs(s$CrI_lower - 110.975) < 1e-10,
            abs(s$CrI_upper - 490.025) < 1e-10,
            s$CrI_upper - s$CrI_lower > 300,
            s$Diagnostic_status == 'review_required', is.na(s$Rhat),
            file.exists(file.path(temp, 'fixture_parameter_diagnostics.csv')),
            file.exists(file.path(temp, 'fixture_acceptance.csv')))
  insufficient <- suppressWarnings(summarize_bedassle_file(path, 'fixture', 200,
    file.path(temp, 'insufficient'), min_saved_states = 50000, make_plots = FALSE))
  stopifnot(!insufficient$summary$Retained_state_target_met)
  make_fixture(path, incomplete = TRUE)
  expect_error(load_bedassle_chain(path, 200))
  make_fixture(path, corrupt = TRUE)
  expect_error(load_bedassle_chain(path, 200))
  expect_error(load_bedassle_chain(file.path(temp, 'missing'), 200))
  cat('BEDASSLE diagnostic and credible-interval tests passed.\n')
}, finally = unlink(temp, recursive = TRUE))
