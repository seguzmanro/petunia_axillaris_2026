# BEDASSLE revision assessment

The original 11 single-chain beta-binomial models used 12,187 LD-filtered SNPs, 2 million iterations and sampling every 250 iterations. Their summary included the full chains and t-based intervals for the mean. Saved-chain checks show substantial autocorrelation and drift in some parameters. The GLMM workflow uses explicit burn-in, coda effective sample sizes and ggmcmc diagnostic reports.

## Authorized changes

1. Retain the same genetic dataset, predictors, priors and proposal scales to assess the effect of a longer run.
2. Set 17,000,000 iterations, 4,250,000 burn-in iterations and sampling every 250 iterations: 68,000 saved states, of which 51,000 are retained for posterior summaries.
3. Calculate 95% equal-tailed posterior credible intervals from retained draws, with numerical endpoints and explicit interval type. Do not use a standard error of the mean as posterior uncertainty.
4. Add coda ESS, autocorrelation, Geweke and Heidelberger-Welch diagnostics, and ggmcmc trace, running-mean, density and autocorrelation reports. Record acceptance rates and flag problematic chains without asserting convergence from iteration count alone.
5. Support re-summarizing old saved runs with explicit burn-in, retaining an audit trail. Write corrected old-run summaries and new long-run outputs under separate prefixes, preserving original files.
6. Validate summary calculations, burn-in boundaries, incomplete-chain handling and the actual model invocation with a small fixture. The user will launch the full workflow on lem_big; do not launch full local runs.

## Risks and execution boundary

Increasing iterations can improve effective sample size but does not improve proposal acceptance or guarantee stationarity. Poorly mixing chains may need tuning and additional runs. A single chain cannot provide between-chain R-hat. Existing non-BEDASSLE working-tree changes must remain untouched. Completion and improved mixing must be checked from the lem_big outputs, not from configuration or launch success.
