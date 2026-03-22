#align(center)[
  #text(size: 20pt, weight: "bold")[Evaluation of Ensemble-Based COBASE -- KNMI Data]

  #v(1em)
  *Jean Nguyen*
]

#v(2em)

We evaluate COBASE methods on real KNMI data with 6 stations, 51 ensemble members, and a 30-day training window. 
$d = 12$ dimensions and $n = 1350$ evaluation days.

== CRPS: EMOS vs EMOS-R

@fig-crps shows the DM test statistics for the CRPS comparing EMOS-Q against EMOS-R across all stations.

#figure(
  image("Results/Figures/flos_dmcrps_t2m_KNMI.pdf"),
  caption: [DM test statistics for CRPS: EMOS-Q vs.\ EMOS-R across KNMI stations.]
) <fig-crps>

== MVPP methods vs COBASE-GCA

@fig-mvpp compares member-based methods (SSh, SimSSh, ECC, GCA) against COBASE-GCA. With the full dataset and 51 ensemble members, ECC can be meaningfully compared to the COBASE variants. ECC performs competitively, achieving an ES of 2.6988 which is close to the best COBASE methods.

#figure(
  image("Results/Figures/flos_dmesvs_mvpp_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: MVPP methods vs.\ COBASE-GCA on KNMI data.]
) <fig-mvpp>


== COBASE shuffling vs parametric copulas

@fig-parcop shows DM test statistics comparing each COBASE variant against its unshuffled parametric copula baseline. COBASE shuffling consistently improves over the raw copula methods for both ES and VS.

#figure(
  image("Results/Figures/flos_dmesvs_parcop_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: COBASE variants vs.\ their parametric copula baselines on KNMI data.]
) <fig-parcop>


== COBASE shuffling improves ensemble copulas

@fig-enscop shows DM test statistics comparing each COBASE-Ens variant against its unshuffled version. As in the synthetic case, shuffling brings clear improvements.

#figure(
  image("Results/Figures/flos_dmesvs_enscop_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: COBASE-Ens variants vs.\ their unshuffled ensemble copula baselines on KNMI data.]
) <fig-enscop>


== Methods vs ECC

@fig-ecc compares member-based and COBASE methods using ECC as the benchmark. Positive values indicate the method is better than ECC.

#figure(
  image("Results/Figures/flos_dmesvs_ecc_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: various methods vs.\ ECC on KNMI data. Positive values indicate the method is better than ECC.]
) <fig-ecc>


== Ensemble-based vs. history-based COBASE

@fig-ensvshist compares ensemble-based COBASE methods against history-based COBASE methods. The history-based variants tend to have slightly better ES, while VS results are more mixed.

#figure(
  image("Results/Figures/flos_dmesvs_ensvshist_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: ensemble-based COBASE vs.\ history-based COBASE on KNMI data. Negative values indicate the history-based variant is better.]
) <fig-ensvshist>


#pagebreak()

== Score summary

@tab-scores reports ES and VS for all methods on the KNMI dataset. The best ES is achieved by COBASE-Frank (2.6905), and the best VS by ECC (298.76). COBASE methods consistently improve over their unshuffled baselines and are competitive with ECC and SimSchaake.

#figure(
  table(
    columns: 3,
    align: (left, right, right),
    stroke: none,
    table.hline(),
    table.header([*Method*], [*ES*], [*VS*]),
    table.hline(),
    [Raw Ensemble],       [3.1251], [354.88],
    table.hline(stroke: 0.3pt),
    [SSh],                [2.7628], [303.34],
    [SimSSh],             [2.6922], [299.17],
    [ECC],                [2.6988], [*298.76*],
    table.hline(stroke: 0.3pt),
    [GCA],                [2.7078], [303.42],
    [Clayton],            [2.7136], [303.21],
    [Frank],              [2.7091], [304.40],
    [Gumbel],             [2.7088], [303.10],
    table.hline(stroke: 0.3pt),
    [COBASE-GCA],         [2.6908], [299.22],
    [COBASE-Clayton],     [2.6920], [299.43],
    [COBASE-Frank],       [*2.6905*], [299.79],
    [COBASE-Gumbel],      [2.6926], [299.34],
    table.hline(stroke: 0.3pt),
    [EnsGCA],             [2.7180], [301.90],
    [EnsClayton],         [2.7261], [301.82],
    [EnsFrank],           [2.7194], [302.32],
    [EnsGumbel],          [2.7216], [302.05],
    table.hline(stroke: 0.3pt),
    [COBASE-EnsGCA],      [2.7019], [298.95],
    [COBASE-EnsClayton],  [2.7069], [299.33],
    [COBASE-EnsFrank],    [2.7022], [299.36],
    [COBASE-EnsGumbel],   [2.7067], [299.51],
    table.hline(),
  ),
  caption: [ES and VS for all methods on KNMI data.],
) <tab-scores>


#pagebreak()

== Computational bottlenecks and numerical issues

Running the full pipeline on the KNMI dataset ($d = 12$, $n = 1350$, $m = 51$) revealed several computational issues that did not appear with the synthetic mock data ($d = 3$, $n approx 300$, $m = 10$).

=== 1. `fitCopula` infinite loops for Archimedean copulas

The `itau` estimator for Archimedean copulas (Frank, Clayton, Gumbel) uses numerical root-finding to invert the Kendall tau-to-parameter relationship. On 3 out of 1350 days, two of the 12 dimensions had *identical rankings* across the 51 ensemble members, producing a pairwise Kendall $tau = 1.0$. The Archimedean copula parameter corresponding to $tau = 1$ is $theta = +infinity$, causing the optimizer to diverge without throwing an error — creating an infinite loop.

*Fix applied:* All `fitCopula` calls are wrapped with `R.utils::withTimeout` (10-second limit). If the fit does not converge, the day falls back to an independence copula.

=== 2. `rMvdc` sampling failures

Even when `fitCopula` succeeds, the resulting copula can produce parameters that cause `rMvdc` (the sampling function) to fail. For the Frank copula specifically, extreme parameters trigger numerical issues in `log1mexp()`, leading to `NaN` values during sampling. This error was not caught by the original `tryCatch` around `fitCopula`.

*Fix applied:* The `rMvdc` call is also wrapped in `tryCatch`, falling back to independent normal samples when sampling fails. Additionally, copula parameters are bounded (Frank: $theta in [0, 50]$, Clayton: $theta in [0, 100]$, Gumbel: $theta in [1, 50]$) to prevent extreme values from reaching the sampler.

=== 3. Parallelization overhead via `future_lapply`

The original implementation used `future_lapply` with `multisession` (16 workers) to parallelize per-day copula fitting. However, each individual `fitCopula` call takes only ~0.02 seconds. The overhead of serializing the large data arrays ($1350 times 51 times 12$) to each worker, combined with `progressr` signaling and worker setup/teardown, far exceeded the computation time. A method that should complete in ~30 seconds sequentially took over 10 minutes with parallel overhead.

*Fix applied:* Replaced `future_lapply` with simple `for` loops, bringing per-method runtime from minutes back to ~30 seconds.

=== Summary of affected days

#figure(
  table(
    columns: 3,
    align: (center, center, left),
    stroke: none,
    table.hline(),
    table.header([*Condition*], [*Days affected*], [*Impact*]),
    table.hline(),
    [$tau = 1.0$ (perfect rank correlation)], [3 / 1350], [Infinite loop in `fitCopula`],
    [$tau >= 0.95$], [6 / 1350], [Slow convergence or failure],
    table.hline(),
  ),
  caption: [Days with degenerate pairwise dependence in the KNMI ensemble.],
) <tab-bottleneck>

These issues arise because Archimedean copulas assume exchangeable dependence (a single parameter governs all pairwise relationships). In 12 dimensions with 2 weather variables per station, this assumption is violated: temperature and dewpoint at the same station can be near-perfectly correlated, while cross-station pairs have weaker dependence.
