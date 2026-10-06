# BEDASSLE runs and posterior diagnostics

Run from this directory on `lem_big`:

```sh
snakemake --snakefile BEDASSLE.snakemake --cores 14
```

Requirements: Snakemake and R packages `vcfR`, `adegenet`, `BEDASSLE`, `foreach`, `doParallel`, `parallel`, `argparse`, `dplyr`, `coda` and `ggmcmc`. The helper `bedassle_diagnostics.R` is tracked as a workflow input. PDF reports use R graphics, with no external TeX compiler.

## Long-run configuration

The existing genetic dataset, 11 champion environmental predictors, default package priors and proposal scales are retained. Each predictor is fitted in a separate beta-binomial model including geographic distance. The configuration specifies:

| Setting | Value |
|---|---:|
| Total iterations | 17,000,000 |
| Burn-in iterations | 4,250,000 (25%) |
| Sampling interval | 250 |
| Total saved states | 68,000 |
| Post-burn-in states | 51,000 |
| Minimum required retained states | 50,000 |
| Checkpoint interval | 100,000 |
| ESS screening threshold | 400 |
| Base random seed | 20261006 |

The interval count and burn-in exclusion are calculated in iteration units. Saved states at or before burn-in are excluded. Configuration validation requires at least 50,000 retained states and a final saved checkpoint. Each predictor receives a distinct reproducible seed (base seed plus predictor index).

The new output prefix is `results/BEDASSLE_champion_noLD_17M`, preserving the original two-million-iteration runs. Existing MCMC output files are never silently overwritten: use a new output prefix for reruns. This runner does not implement checkpoint continuation. Intermediate checkpoints are retained by BEDASSLE but cannot be summarized as completed chains.

## Outputs and interpretation

- `_BEDASSLE_RES_CI.csv`: mean and median of retained αE/αD draws, numerical lower and upper bounds of a 95% equal-tailed posterior credible interval (2.5th and 97.5th percentiles), retained-state count, ESS, and diagnostic flags. The legacy `CI` column contains a formatted credible interval, explicitly identified by `CI_type`; the numerical bounds retain full precision.
- `_parameter_diagnostics.csv`: per-predictor summaries and diagnostics for α0, αD, αE, α2, β, population-specific φ and αE/αD.
- `_mcmc_plots/`: per-model ESS tables, proposal acceptance rates, numeric diagnostic tables, ggmcmc posterior density/trace/running-mean/autocorrelation/Geweke reports, and full-chain traces showing the burn-in boundary. The ggmcmc input explicitly disables additional automatic burn-in removal.
- Per-predictor `_chain.log`: iteration progress, model seed and elapsed runtime.
- `_run_metadata.rds` and `_session_info.txt`: input checksums, settings, package versions and diagnostic scope. Summary-only mode records its settings separately.

`review_required` flags ESS below 400, unavailable diagnostics, |Geweke Z| > 1.96, a failed Heidelberger-Welch stationarity/halfwidth check, suggested additional within-test discard, or insufficient retained states. These are screening flags for human review; individual tests can disagree and repeated testing can flag otherwise satisfactory parameters. `no_numeric_flags` is not proof of convergence: inspect trace, running-mean and autocorrelation plots before interpretation. No additional states are silently discarded in response to diagnostic tests.

The original single-chain-per-model design is retained to assess the impact of longer runs. R-hat is reported as unavailable because it requires multiple chains for the same model; chains from different predictors must not be combined. BEDASSLE does not save locus-specific θ and μ trajectories, so those trajectories cannot be diagnosed retrospectively; proposal acceptance summaries for those families are provided. Longer runs do not improve proposal acceptance by themselves. If chains still drift or ESS remains low, tuning and independent replicate chains should be considered before using results in the manuscript.

## Correct the summaries of existing runs without rerunning MCMC

For the original two-million-iteration outputs, explicitly specify an appropriate exploratory burn-in and a new output prefix:

```sh
Rscript bedassle_script.R \
  --summarize_only \
  --summary_input_prefix results/BEDASSLE_champion_noLD \
  --env_dist ../EnvironValriables/champion_EnvironVars_Selected_Distances.csv \
  --out_prefix results/BEDASSLE_champion_noLD_corrected_legacy \
  --burnin 500000
```

This leaves 6,000 retained states. It generates correctly defined posterior interval calculations but flags the 50,000-state target as unmet. The 500,000-iteration exclusion is a diagnostic choice, not evidence that the original chains converged. The resulting intervals remain provisional if mixing is poor. Use `--no_plots` for numerical-only post-processing. Missing, corrupt or incomplete chains cause an error rather than a silently incomplete final table.

After a completed long run, the same mode can regenerate reports using the new chain prefix and a distinct summary output prefix; set `--burnin 4250000` for the configured long runs.

## Verification

```sh
Rscript tests/test_bedassle_diagnostics.R
snakemake --snakefile BEDASSLE.snakemake --cores 14 --dry-run
```

The deterministic R checks cover iteration accounting, exact burn-in boundaries, posterior distribution intervals, retained-state flags and rejection of missing/corrupt/incomplete outputs. Small local MCMC fixtures can validate execution/report generation; they do not validate biological inference or mixing of the full-data long runs.

Validation for this revision: deterministic R checks passed; Snakemake 9.27.0 dry-run resolved the configured inputs and both aggregate outputs; an eight-locus, two-predictor fixture completed through the actual parallel MCMC runner; summary-only mode generated numeric diagnostics and PDF reports; all 11 original full-data saved chains were reprocessed under a separate corrected-legacy prefix. Full-data 17-million-iteration runs are left for execution on `lem_big`.
