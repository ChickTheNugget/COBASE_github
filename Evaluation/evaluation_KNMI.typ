#align(center)[
  #text(size: 20pt, weight: "bold")[Evaluation of Ensemble-Based COBASE -- KNMI Data]

  #v(1em)
  *Jean Nguyen*
]

#v(2em)

We evaluate COBASE methods on real KNMI data with 6 stations, 51 ensemble members, and a 30-day training window.
$d = 12$ dimensions and $n = 1350$ evaluation days.

/*
== Univariate postprocessing: EMOS-Q vs EMOS-R

@fig-crps shows DM test statistics for the CRPS, comparing EMOS-Q against the EMOS-R benchmark per station for both temperature (T2m) and dew point temperature (DPT).

#figure(
  image("../Results/Figures/flos_dmcrps_t2m_KNMI.pdf"),
  caption: [DM test statistics for CRPS: EMOS-Q vs.\ EMOS-R per station for T2m and DPT on KNMI data.]
) <fig-crps>

== SimSSh vs SimSSh-R

@fig-simssh shows DM test statistics for the Energy Score and Variogram Score, comparing SimSSh against SimSSh-R.

#figure(
  image("../Results/Figures/flos_dmesvs_simssh_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: SimSSh vs.\ SimSSh-R on KNMI data.]
) <fig-simssh>

== MVPP methods vs COBASE-GCA

@fig-mvpp compares standard multivariate postprocessing methods against COBASE-GCA as benchmark.

#figure(
  image("../Results/Figures/flos_dmesvs_mvpp_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: MVPP methods vs.\ COBASE-GCA on KNMI data.]
) <fig-mvpp>

== COBASE shuffling improves parametric copulas

@fig-parcop shows DM test statistics comparing each COBASE variant against its unshuffled parametric copula baseline.

#figure(
  image("../Results/Figures/flos_dmesvs_parcop_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: COBASE variants vs.\ their unshuffled parametric copula baselines on KNMI data.]
) <fig-parcop>

== COBASE shuffling improves ensemble copulas

@fig-enscop shows DM test statistics comparing each COBASE-Ens variant against its unshuffled version.

#figure(
  image("../Results/Figures/flos_dmesvs_enscop_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: COBASE-Ens variants vs.\ their unshuffled ensemble copula baselines on KNMI data.]
) <fig-enscop>

As in the synthetic case, shuffling brings clear improvements.

== Methods vs ECC

@fig-ecc compares history-based methods against ECC. Positive values indicate the method outperforms ECC.

#figure(
  image("../Results/Figures/flos_dmesvs_ecc_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: history-based methods vs.\ ECC on KNMI data. Positive values indicate the method outperforms ECC.]
) <fig-ecc>
*/

== Ensemble-based COBASE vs ECC

@fig-ecc-ens compares the ensemble-based COBASE methods against ECC. Positive values indicate the method outperforms ECC.

#figure(
  image("../Results/Figures/flos_dmesvs_ecc_ens_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: ensemble-based COBASE methods vs.\ ECC on KNMI data. Positive values indicate the method outperforms ECC.]
) <fig-ecc-ens>

== All methods vs ECC

@fig-ecc-combined shows all COBASE methods (history-based and ensemble-based) compared against ECC in a single view.

#figure(
  image("../Results/Figures/flos_dmesvs_ecc_combined_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: all methods vs.\ ECC on KNMI data. Positive values indicate the method outperforms ECC.]
) <fig-ecc-combined>


== Ensemble-based vs. history-based COBASE

@fig-ensvshist compares ensemble-based COBASE methods against history-based COBASE methods. The history-based variants tend to have slightly better ES, while VS results are more mixed.

#figure(
  image("../Results/Figures/flos_dmesvs_ensvshist_KNMI.pdf"),
  caption: [DM test statistics for ES and VS: ensemble-based COBASE vs.\ history-based COBASE on KNMI data. Negative values indicate the history-based variant is better.]
) <fig-ensvshist>


#pagebreak()

== Score summary

@tab-scores reports ES and VS for all methods on the KNMI dataset.

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
