#align(center)[
  #text(size: 20pt, weight: "bold")[Evaluation of Ensemble-Based COBASE]
  
  #v(1em)
  *Jean Nguyen*
]

#v(2em)

We fit copulas on the current day's ensemble members rather than on historical observations. Experiments are run on mock data with 3 stations, 10 ensemble members, and a 30-day training window.

== Copula fitting difference


*1. History-based Fitting:* The copula is trained on a rolling window of the past 30 days of historical observations.
#figure(
  table(
    columns: 4,
    align: center,
    stroke: none,
    table.hline(),
    table.header([*Date*], [*St 1*], [*St 2*], [*St 3*]),
    table.hline(),
    [2021-01-02], [276.74], [275.84], [277.16],
    [2021-01-03], [276.36], [273.89], [274.62],
    [2021-01-04], [273.03], [274.23], [272.78],
    [...], [...], [...], [...],
    [2021-01-29], [280.29], [280.27], [276.14],
    [2021-01-30], [271.41], [270.64], [269.45],
    [2021-01-31], [269.90], [270.17], [271.25],
    table.hline(),
  ),
  caption: [History-based training data: 30 past historical observations],
) <tab-hist-data>

*2. Ensemble-based Fitting:* The copula is trained on the current day's raw ensemble members.
#figure(
  table(
    columns: 11,
    align: center,
    stroke: none,
    table.hline(),
    table.header([*Station*], [*$M_0$*], [*$M_1$*], [*$M_2$*], [*$M_3$*], [*$M_4$*], [*$M_5$*], [*$M_6$*], [*$M_7$*], [*$M_8$*], [*$M_9$*]),
    table.hline(),
    [1], [270.711], [270.860], [272.666], [269.897], [270.639], [273.253], [271.748], [273.485], [271.426], [271.420],
    [2], [270.918], [273.859], [270.974], [272.195], [272.360], [272.021], [274.175], [273.721], [274.034], [274.171],
    [3], [266.830], [268.034], [268.870], [269.126], [266.778], [268.380], [268.310], [269.012], [265.745], [268.502],
    table.hline(),
  ),
  caption: [Ensemble-based training data: 10 current ensemble members]
) <tab-ens-data>


/*
== Univariate postprocessing: EMOS-Q vs EMOS-R

@fig-crps shows DM test statistics for the CRPS, comparing EMOS-Q against the EMOS-R benchmark per station for both temperature (T2m) and dew point temperature (DPT).

#figure(
  image("../Results/Figures/flos_dmcrps_t2m_Mock_data.pdf"),
  caption: [DM test statistics for CRPS: EMOS-Q vs.\ EMOS-R per station for T2m and DPT on mock data.]
) <fig-crps>

== SimSSh vs SimSSh-R

@fig-simssh shows DM test statistics for the Energy Score and Variogram Score, comparing SimSSh against SimSSh-R.

#figure(
  image("../Results/Figures/flos_dmesvs_simssh_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: SimSSh vs.\ SimSSh-R on mock data.]
) <fig-simssh>

== MVPP methods vs COBASE-GCA

@fig-mvpp compares standard multivariate postprocessing methods against COBASE-GCA as benchmark.

#figure(
  image("../Results/Figures/flos_dmesvs_mvpp_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: MVPP methods vs.\ COBASE-GCA on mock data.]
) <fig-mvpp>

== COBASE shuffling improves parametric copulas

@fig-parcop shows DM test statistics comparing each COBASE variant against its unshuffled parametric copula baseline.

#figure(
  image("../Results/Figures/flos_dmesvs_parcop_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: COBASE variants vs.\ their unshuffled parametric copula baselines on mock data.]
) <fig-parcop>

== COBASE shuffling improves ensemble copulas

@fig-enscop shows DM test statistics comparing each COBASE-Ens variant against its unshuffled version.

#figure(
  image("../Results/Figures/flos_dmesvs_enscop_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: COBASE-Ens variants vs.\ their unshuffled ensemble copula baselines.]
) <fig-enscop>

The @fig-enscop points at shuffling being beneficial, showing consistent improvements over non shuffled baselines.

== Methods vs ECC

@fig-ecc compares history-based methods against ECC. Positive values indicate the method outperforms ECC.

#figure(
  image("../Results/Figures/flos_dmesvs_ecc_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: history-based methods vs.\ ECC on mock data. Positive values indicate the method outperforms ECC.]
) <fig-ecc>
*/

== Ensemble-based COBASE vs ECC

@fig-ecc-ens compares the ensemble-based COBASE methods against ECC. Positive values indicate the method outperforms ECC.

#figure(
  image("../Results/Figures/flos_dmesvs_ecc_ens_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: ensemble-based COBASE methods vs.\ ECC on mock data. Positive values indicate the method outperforms ECC.]
) <fig-ecc-ens>

== All methods vs ECC

@fig-ecc-combined shows all COBASE methods (history-based and ensemble-based) compared against ECC in a single view.

#figure(
  image("../Results/Figures/flos_dmesvs_ecc_combined_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: all methods vs.\ ECC on mock data. Positive values indicate the method outperforms ECC.]
) <fig-ecc-combined>


== Ensemble-based vs. history-based COBASE

@fig-ensvshist compares ensemble-based COBASE methods against history-based COBASE methods. For the Energy Score, the history-based variants are significantly better across all copula families. For the Variogram Score, results are mixed: COBASE-EnsFrank and COBASE-EnsGumbel are comparable to their history-based counterparts, while COBASE-EnsGCA performs worse.

#figure(
  image("../Results/Figures/flos_dmesvs_ensvshist_Mock_data.pdf"),
  caption: [DM test statistics for ES and VS: ensemble-based COBASE vs.\ history-based COBASE. Negative values indicate the history-based variant is better.]
) <fig-ensvshist>


#pagebreak()

== Score summary

@tab-scores reports ES and VS for all methods. 

#figure(
  table(
    columns: 3,
    align: (left, right, right),
    stroke: none,
    table.hline(),
    table.header([*Method*], [*ES*], [*VS*]),
    table.hline(),
    [Raw Ensemble],       [1.5686], [9.4873],
    table.hline(stroke: 0.3pt),
    [SimSSh],             [1.5924], [9.3181],
    [ECC],                [1.5916], [9.1604],
    table.hline(stroke: 0.3pt),
    [GCA],                [1.6621], [10.1238],
    [Clayton],            [1.6746], [10.1746],
    [Frank],              [1.6340], [9.8692],
    [Gumbel],             [1.6573], [9.7173],
    table.hline(stroke: 0.3pt),
    [COBASE-GCA],         [1.5941], [9.2952],
    [COBASE-Clayton],     [1.5929], [9.2506],
    [COBASE-Frank],       [1.5907], [9.2940],
    [COBASE-Gumbel],      [*1.5902*], [9.3514],
    table.hline(stroke: 0.3pt),
    [EnsGCA],             [1.6654], [10.5961],
    [EnsClayton],         [1.6776], [10.3395],
    [EnsFrank],           [1.6731], [9.9944],
    [EnsGumbel],          [1.6674], [9.8005],
    table.hline(stroke: 0.3pt),
    [COBASE-EnsGCA],      [1.6021], [9.3977],
    [COBASE-EnsClayton],  [1.5979], [9.2343],
    [COBASE-EnsFrank],    [1.5946], [*9.1810*],
    [COBASE-EnsGumbel],   [1.5928], [9.2667],
    table.hline(),
  ),
  caption: [ES and VS for all methods on mock data.],
) <tab-scores>


== Algorithm Implementation Summary

The core algorithmic difference implemented in the codebase lies purely in the data subset used to estimate the copula parameters. Below is a simplified representation of the multivariate postprocessing workflow for a given target day:

```text
For each target forecast day t:

  // 1. Univariate Postprocessing
  U_ens = apply_EMOS(raw_ensemble[t])  // Correct biases per station

  // 2. Copula Parameter Estimation
  IF method == "History-based":
      train_data = observations[t-30 to t-1]
  ELSE IF method == "Ensemble-based":
      train_data = raw_ensemble_members[t]
  
  theta = estimate_copula_parameters(train_data)
  fitted_copula = define_copula(family, theta)

  // 3. Rank Template Generation
  copula_samples = sample_from(fitted_copula, N_members)
  rank_template = extract_ranks(copula_samples)

  // 4. COBASE Shuffling
  final_ensemble = reorder_array(U_ens, rank_template)
```
