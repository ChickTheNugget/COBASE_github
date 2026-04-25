#align(center)[
  #text(size: 20pt, weight: "bold")[Comparing Copula Parameter Fits: Ensemble vs. Observation]

  #v(1em)
  *Jean Nguyen*
]

#v(2em)


We fit an Archimedean copula (Clayton, Frank, or Gumbel) in two different ways:

- *Obs-based*: the copula is fit on the past 30 days of observations
- *Ens-based*: the copula is fit on the current day's ensemble members (10 for Mock_data, 51 for KNMI)

Each fit produces a single dependence parameter $theta$ for that day. 
For each day, we now have 2 $theta$ values that, are estimating the same quantity (the dependence among stations on that day). The four analyses that follow ask: do they agree, by how much do they differ, and what does each one's behaviour look like overall?


Just like we used Kendall's tau to estimate the Copulas we can use the parameter values to get the Kendall's tau coefficients.

#align(center)[
  $tau_("Clayton") = theta / (theta + 2), quad tau_("Gumbel") = 1 - 1/theta, quad tau_("Frank")$ via the Debye function.
]

== 1. Scatter plot: $theta_("ens")$ vs $theta_("obs")$

#figure(
  image("../Results/Figures/ParamComparison/scatter_Mock_data.pdf", width: 100%),
  caption: [Scatter of $theta_("ens")$ vs $theta_("obs")$ on Mock_data, per copula family.]
) <fig-scatter-mock>

#figure(
  image("../Results/Figures/ParamComparison/scatter_KNMI.pdf", width: 100%),
  caption: [Scatter of $theta_("ens")$ vs $theta_("obs")$ on KNMI, per copula family.]
) <fig-scatter-knmi>

- Points are consistently above the line  $arrow.r$ ens-based $theta$ is systematically larger.
- Mock_data: $|r| lt.eq 0.04$, $|rho| lt.eq 0.13$ $arrow.r$ essentially uncorrelated day-to-day.
- KNMI: $r approx 0.16$–$0.19$, $rho approx 0.17$–$0.20$ $arrow.r$ weakly aligned, far from one-to-one.

== 2. Bland-Altman plot

#figure(
  image("../Results/Figures/ParamComparison/blandaltman_Mock_data.pdf", width: 100%),
  caption: [Bland-Altman plot on Mock_data, per copula family.]
) <fig-ba-mock>

#figure(
  image("../Results/Figures/ParamComparison/blandaltman_KNMI.pdf", width: 100%),
  caption: [Bland-Altman plot on KNMI, per copula family.]
) <fig-ba-knmi>


- Bias line is positive in every panel $arrow.r$ confirms ens overestimates $theta$.
- Bias is small for Gumbel (lower-bounded by 1), much larger for Clayton and Frank.
- 95% bands are wide, especially Frank $arrow.r$ daily disagreement of several units of $theta$ is common.
- KNMI: clear diagonal trend $arrow.r$ proportional bias (ens $approx$ scaled-up obs, not just shifted).

== 3. Distribution of $theta$ over the test period

#figure(
  image("../Results/Figures/ParamComparison/distribution_Mock_data.pdf", width: 100%),
  caption: [Distribution of $theta$ over the test period on Mock_data.]
) <fig-dist-mock>

#figure(
  image("../Results/Figures/ParamComparison/distribution_KNMI.pdf", width: 100%),
  caption: [Distribution of $theta$ over the test period on KNMI.]
) <fig-dist-knmi>


- Obs-based boxes are narrow; ens-based boxes are much wider with more outliers.
- Reason: 30-day window smooths the obs-based estimate; ens-based uses only one day of members so it's noisier.
- Ens-based medians sit higher than obs-based $arrow.r$ confirms the upward bias.

== 4. Summary table with implied Kendall's $tau$

#figure(
  table(
    columns: 7,
    align: (left, right, right, right, right, right, right),
    stroke: none,
    table.hline(),
    table.header(
      [*Family*], [*n*],
      [*mean $theta_("obs")$*], [*mean $theta_("ens")$*],
      [*mean $tau_("obs")$*], [*mean $tau_("ens")$*],
      [*mean diff $theta$*]
    ),
    table.hline(),
    [Clayton], [240], [0.296], [0.486], [0.124], [0.174], [0.201],
    [Frank],   [240], [1.101], [1.609], [0.119], [0.166], [0.514],
    [Gumbel],  [240], [1.153], [1.216], [0.128], [0.158], [0.064],
    table.hline(),
  ),
  caption: [Summary statistics on Mock_data (10 ensemble members, 240 forecast days).]
) <tab-summary-mock>

#figure(
  table(
    columns: 7,
    align: (left, right, right, right, right, right, right),
    stroke: none,
    table.hline(),
    table.header(
      [*Family*], [*n*],
      [*mean $theta_("obs")$*], [*mean $theta_("ens")$*],
      [*mean $tau_("obs")$*], [*mean $tau_("ens")$*],
      [*mean diff $theta$*]
    ),
    table.hline(),
    [Clayton], [1350], [0.250], [0.773], [0.108], [0.252], [0.523],
    [Frank],   [1350], [0.749], [2.268], [0.082], [0.226], [1.549],
    [Gumbel],  [1350], [1.161], [1.404], [0.136], [0.264], [0.244],
    table.hline(),
  ),
  caption: [Summary statistics on KNMI (51 ensemble members, 1350 forecast days).]
) <tab-summary-knmi>

- *Bias.* Ens-based fitting always returns a larger mean $theta$ than obs-based fitting. The bias on Kendall's $tau$ is roughly $+ 0.05$ on Mock_data and roughly $+ 0.13$ on KNMI.
- *Cross-family consistency.* All three families tell the same story (positive bias, more noise in ens).

== Conclusion

The two methods do not give the same $theta$. 

1. *They rarely agree on the same day.* On Mock_data the paired estimates are basically uncorrelated; on KNMI there is a small positive correlation ($r approx 0.17$). One estimate tells you very little about the other.
2. *Ens-based always gives a larger $theta$.* In every family on every dataset, ens-based estimates more dependence than obs-based. The gap on Kendall's $tau$ ranges from about $+0.03$ (Mock Gumbel) to $+0.14$ (KNMI Clayton/Gumbel).
3. *Ens-based is much noisier.* Wider boxplots, more outliers, and more failed fits — especially when the ensemble is small.
