# Statistical calibration and secondary analyses

The primary result is theta = 0.42096778795400724, one-sided competitive p = 0.7438592561407439. The calibration and secondary calculations were specified in `protocol/revision_analysis_plan.json` after the primary result and preliminary secondary estimates were known, and before the specified outputs were generated. Their interpretation is secondary or exploratory. All results are computed from retained scientific inputs.

## 1. Estimand and uncertainty

Theta is a normalized Mann-Whitney U statistic, with half weight for ties. In this finite table it is the fraction of candidate-background gene pairs favoring the candidate: approximately 42.1%. It is not the probability that a modification protects a cell. The original one-sided random-set reference samples six distinct genes uniformly from the same 11,689 reported point estimates.

The observed rank shift is opposite to the hypothesized direction, and the competitive test did not detect positive enrichment. Neither significant depletion nor equivalence was tested or established. Gene-label simulation quantifies a reference distribution conditional on this table. It does not measure uncertainty in the source effect estimates, their ordering under repeated experiments, or transfer to another biological model. The available exported scores and adjusted q-values do not supply the original nominal p-values or joint effect covariance. No nominal p-values are reconstructed from q-values, and no six-gene bootstrap is labeled biological uncertainty.

The normalized U identification follows [Mann and Whitney (1947)](https://doi.org/10.1214/aoms/1177730491). The separation of a nonsignificant result from evidence of absence, and of simulation precision from scientific uncertainty, follows the interpretation cautions in [Greenland et al. (2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4877414/).

## 2. Original six-gene null calibration

We repeated the original PCG64 seed 20261025 and 1,000,000 draws using the identical batch size and whole-row collision rejection. The primary exceedance count and p-value reproduce exactly. Empirical quantiles use the `higher` rule, which selects an observed simulated statistic rather than interpolating between values.

| Reference quantity | Independently recomputed value |
|---|---:|
| Observed theta | 0.420968 |
| Null mean | 0.500014 |
| Null standard deviation | 0.117956 |
| Central 95% lower boundary | 0.270778 |
| Central 95% upper boundary | 0.730035 |
| Empirical upper 5% boundary | 0.694827 |
| Fraction of null draws strictly above that boundary | 0.049993 |
| Draws exactly at that boundary | 14 |

These are competitive-null reference quantities, **not confidence limits for biological effects**. The observed value lies comfortably within this broad null distribution. The reference rejection rule used in the design-sensitivity simulations is theta strictly greater than 0.694827241861394. The boundary is estimated from the retained Monte Carlo sample; it is not claimed to be an exact enumerated finite-universe critical value. Strict inequality makes boundary ties explicit.

One million draws precisely estimate the reference tail probability. They do not change the fact that only six genes were nominated. The entire draw vector, fixed-bin histogram and empirical quantiles are retained for plotting and independent checking.

## 3. Why zero source-threshold hits is not a second strong negative result

Exactly 98 of the 11,689 genes have both a positive signed score and reported genome-wide q <= 0.05. Under uniform six-gene sampling without replacement,

`P(no source-threshold hits) = choose(11689 - 98, 6) / choose(11689, 6) = 0.9507286242496122`.

The expected count is only 0.050304. Zero observed hits is therefore common under the same unrestricted random-set reference. It remains correct to report that none of the six meets this source threshold, but this is not independent strong evidence that their biological effects are absent. The source's genome-wide discovery q-values are retained. A justified targeted analysis would require original nominal statistics or a suitable experimental model; adjusted q-values cannot be inverted to provide those missing inputs.

## 4. Secondary joint conditional matching

### Question and rules

This analysis asks whether the six nominees outrank comparable genes under the already declared guide-coverage and abundance rules. It uses the 11,474 scored genes with mapped raw covariates; the other 215 remain explicitly outside this reference. All six nominees are retained. Strata are exact triples of intended guide count, no-PBMC mean-CPM quintile, and pretransplant mean-CPM quintile. The original stored quintiles were independently checked against the original `pandas.qcut` rule.

Let stratum s contain n_s genes and k_s nominees. Let theta_s compare those k_s nominees with the remaining n_s-k_s genes, using half weight for ties. The matched statistic is

`T = sum_s (k_s / 6) * theta_s`.

The conditional null preserves k_s and uniformly samples that many distinct genes from the **entire stratum, including original nominees**. Exact strata partition the universe, so genes cannot collide across strata. Sampling within a stratum is without replacement; whole colliding tuples are rejected, preserving the intended uniform subset law. Finite-population dependence within each stratum is retained. Biological correlation among gene effects is neither estimated nor corrected by this procedure.

In this dataset the six nominees occupy six different strata and each k_s = 1. Thus observed T equals the mean of the six original descriptive matched percentiles. The original reported peer counts excluded nominees; the conditional sampling pools include them.

| Nominee | Full stratum size | Noncandidate peers | Within-stratum theta |
|---|---:|---:|---:|
| RNLS | 1,696 | 1,695 | 0.336873 |
| HIVEP2 | 543 | 542 | 0.647601 |
| ZNF28 | 1,215 | 1,214 | 0.126853 |
| MED11 | 557 | 556 | 0.097122 |
| ZBED3 | 434 | 433 | 0.062356 |
| GCN1 | 347 | 346 | 0.942197 |

### Result and limits

Observed T = **0.36883375029516285**. In 1,000,000 PCG64 draws with seed 20261026, 863,280 null statistics were at least as large, giving plus-one one-sided p = **0.8632801367198633** and Monte Carlo SE = **0.0003435514**. No positive conditional enrichment was detected. This is a different, within-stratum estimand and reference from the unrestricted primary theta; the numerical change from 0.421 to 0.369 is not an adjusted causal effect.

GCN1's higher matched position and HIVEP2's more moderate one remain candidate-specific descriptive distinctions. They do not make GCN1 a validated edit or provide a newly computed targeted gene p-value. Matching measured abundance and guide coverage does not establish exchangeability or comparable editing. Control graft abundance is post-transplant and can itself reflect biological selection. Gene correlations remain a separate concern: [Wu and Smyth's CAMERA study (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/) explains why competitive tests can be affected by them. CAMERA was **not** applied to the published vector; the appropriate experimental inputs and covariance model for that analysis are not available in this export.

## 5. Common-universe model comparison and reference centering

Published-score and raw-count-derived estimates were compared using the same 11,474 mapped genes and the same six nominees. No separate score-specific universe was selected.

| Score definition | Theta | Positive nominee scores |
|---|---:|---:|
| Published signed score | 0.42086675967910708 | 2 |
| Descriptive log2(CPM + 1) contrast | 0.49504418090919661 | 0 |
| Same descriptive contrast minus shared NT contrast | 0.49504418090919661 | 2 |

The common non-targeting contrast is -0.5442449920913339. If `d'_g = d_g - c`, then `d'_g > d'_h` exactly when `d_g > d_h`, and ties are also preserved. Therefore every pairwise contribution and theta are unchanged when the universe is fixed. The complete computed rank vectors are identical before and after centering, not merely the six candidate ranks.

Restored sign agreement after non-targeting centering explains the reference zero. It is not additional confirmation of the primary rank-enrichment result. The near-neutral raw theta and lower published-score theta instead show dependence on the score definition. The descriptive CPM contrast does not reproduce the source's negative-binomial, non-targeting-referenced mixed-model interaction, so this is an exploratory estimand comparison rather than an independent biological replication or a refutation of that model.

## 6. Sensitivity under declared stylized rank alternatives

These simulations answer which *specified rank patterns* this six-gene design would tend to detect. Their parameters were not estimated from the observed theta. They are not observed power, inferred biological effects, or plausible-effect claims.

For nominee position i, independently draw a rank quantile `U_i ~ Beta(alpha_i, 1)`, then assign zero-based finite-universe rank `j_i = floor(11689 * U_i)`. Before distinctness conditioning, the probability of selecting rank j is

`p_ij = ((j + 1) / 11689)^alpha_i - (j / 11689)^alpha_i`.

If any two ranks coincide, the entire six-rank tuple is resampled. This yields the independent product rank law **conditional on distinct genes**. Diffuse alternatives give all six positions the same alpha. Sparse alternatives give one or three positions alpha = 5 and the remaining positions alpha = 1. Positions label sampling roles, not particular genes. The symmetric statistic makes which position carries a role immaterial. Higher alpha increases preference for higher source ranks; no mechanism or expression model is implied.

Each scenario uses 200,000 sets and its frozen PCG64 seed. Detection means theta exceeds the empirical primary-null boundary 0.694827. The nominal quantile mean below is the average `alpha_i/(alpha_i+1)` before finite-rank discretization and distinctness conditioning; the actual simulated theta mean is reported separately. These means must not be confused with the measured biological effect of an edit.

| Stylized scenario | Nominal mean rank quantile | Simulated mean theta | Rejection fraction | MC SE of fraction |
|---|---:|---:|---:|---:|
| Uniform null, all alpha = 1 | 0.500 | 0.500641 | 0.049955 | 0.000487 |
| Diffuse, all alpha = 1.5 | 0.600 | 0.600358 | 0.194400 | 0.000885 |
| Diffuse, all alpha = 2 | 0.667 | 0.666651 | 0.400855 | 0.001096 |
| Diffuse, all alpha = 3 | 0.750 | 0.749800 | 0.763190 | 0.000951 |
| Diffuse, all alpha = 5 | 0.833 | 0.833580 | 0.984165 | 0.000279 |
| One alpha = 5, five alpha = 1 | 0.556 | 0.555260 | 0.107040 | 0.000691 |
| Three alpha = 5, three alpha = 1 | 0.667 | 0.666804 | 0.388215 | 0.001090 |

The uniform scenario's rejection fraction agrees with the nominal reference level. Moderate diffuse rank preference or one strongly preferred rank often does not cross the six-gene threshold in these stylized models. A sufficiently strong diffuse pattern is much easier to detect. Similar mean theta does not uniquely determine detection probability: the six alpha = 2 scenario and the three alpha = 5 scenario have different rank distributions and different rejection fractions. Unmodeled biological correlation, effect-estimation error and another alternative distribution could change these results. Monte Carlo SE quantifies simulation noise conditional on the frozen threshold and scenario, not uncertainty that these scenarios describe biology.

## Reproduction and artifact dictionary

`python src/revision_analysis.py` verifies strict frozen protocol hashes and the content of regenerated tables against the frozen revision references and rebuilds all `results/revision_*` analysis files. It imports the original rank/sampling functions without modifying `src/transfer.py`. `python -m unittest tests.test_revision_analysis -v` checks small-universe enumeration, tie handling, exact stratified-statistic calculations, the weighted collision-rejection law, hypergeometric identities, constant-shift invariance and unchanged primary results.

* `revision_summary.json`: exact calibration, no-hit probability, matched test and invariance findings.
* `revision_primary_null_histogram.tsv` and `revision_primary_null_quantiles.tsv`: plotting data with clear reference labels.
* `revision_primary_null.npy.gz` and `revision_matched_null.npy.gz`: deterministic compressed complete draw vectors (load with NumPy after gzip decompression).
* `revision_matched_strata.tsv` and `revision_matched_null_histogram.tsv`: conditional-reference membership accounting and distribution.
* `revision_common_universe.tsv`: paired score-definition comparison on one denominator.
* `revision_alternative_sensitivity.tsv`: alpha vectors, seeds, simulated theta means, rejection counts/fractions, MC errors and collision counts.

Set-level nondetection, zero discovery-threshold hits and absent post-immune function measurements are distinct observations, not three independent biological-failure findings.
