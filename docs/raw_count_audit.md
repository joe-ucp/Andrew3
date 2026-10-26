# Raw-count robustness and denominator audit

The raw-data arm supports a bounded interpretation of the published-table result. It does not reconstruct the paper's fitted model or create independent biological validation. The public count file has 77,441 unique guide barcodes, matches the retained source SHA-256 exactly, and joins one-to-one to the GEO construct reference and Broad original-target gene annotation. All six primary candidates have four intended guides. Every graft sample remains in analysis.

## Source and sample identity

`evidence/raw_count_sources.json` freezes URLs, retrieval times, source bytes and SHA-256. `src/fetch_raw_counts.py` retrieves the same public sources and rejects changed bytes against the frozen manifest. The main guide annotation is Broad CP0043 original-target chip dated 2022-03-24; current-sequence annotation is independently audited against Broad GRCh38 NCBI strict chip dated 2024-01-30. The Addgene Brunello contents are also retained. These sources establish guide identity; GEO's CP0043 reference alone maps sequences to constructs, not genes.

GEO GSM6008413 sample information explicitly identifies no-PBMC graft mice 64, 69, 70, 73, 74 and 78, and hPi graft mice 65, 67, 68, 72, 76 and 79. `DA446 Pre-Tx` and `CP0043 Library DNA` are baseline materials, not independent graft mice. `NTC` and the unlabeled zero column are excluded from normalization and contrasts. No donor labels are supplied, so neither six independent donors nor donor-adjusted effects are claimed.

Sequencing coverage is uneven: control mouse 69 has 2,499,075 reads and 45.83% zero guides; hPi mouse 79 has only 42,430 reads and 77.89% zero guides. Other graft samples contain 4.46–17.07 million reads. No sample was removed after inspecting outcomes. The full QC table and leave-one-animal-out contrasts retain this limitation.

## Published score semantics

Sintov supplementary `mmc1.pdf`, “Whole genome CRISPR screen in vivo and analysis” (model on PDF page 28), models guide counts using `glmer.nb` from lme4, with fixed KO, PBMC, and KO-by-PBMC terms and random intercepts for mouse and guide. The stated model is `log(counts) ~ 1 + KO + PBMC + KO*PBMC + (1|Mouse) + (1|Guide)`. Interaction significance uses a likelihood-ratio comparison with an additive model; p-values are adjusted using Benjamini–Hochberg FDR. The published gene universe excludes genes not expressed in at least 10% of either SC-alpha or SC-beta cells in any assessed condition. Main Figure 1G describes enrichment relative to non-targeting guides. [Primary article and supplementary information](https://pmc.ncbi.nlm.nih.gov/articles/PMC9481918/).

These materials establish a signed KO-by-environment contrast and its protective direction. Main Figure 2E's caption explicitly describes positively enriched CRISPR screen hits using a `log2 fold change >1` threshold, in addition to transcriptomic thresholds. That label suggests a log2 screen scale. However, the released workbook and model description do not explicitly connect the workbook field named `Fold change` to a particular exported coefficient transformation. Therefore this project preserves the workbook values and uses ranks, without exponentiating the field or calling its magnitude a conventional fold ratio. The descriptive raw contrasts below are not numerical estimates of that fitted interaction.

## Prespecified descriptive analysis

`protocol/raw_sensitivity.json` was written before candidate raw outcomes were inspected. For each sample, counts are divided by total counts and multiplied by one million; each guide becomes `log2(CPM+1)`. Guides are averaged within each gene and mouse, then the six hPi mouse means are compared with the six no-PBMC mouse means. Mouse omissions recompute means with the other 11 animals; guide omissions recompute the within-mouse gene summaries. Omission ranges are stability diagnostics, not confidence intervals or independent replicate counts.

| Gene | Mean log2(CPM+1) contrast | Animal omission range | Guide omission range |
|---|---:|---:|---:|
| HIVEP2 | -0.104 | -0.633 to 0.174 | -0.345 to 0.202 |
| RNLS | -0.648 | -1.090 to -0.336 | -0.897 to -0.335 |
| GCN1 | -0.211 | -0.446 to 0.489 | -0.337 to -0.033 |
| MED11 | -0.792 | -1.075 to -0.401 | -1.030 to -0.543 |
| ZBED3 | -0.880 | -1.103 to -0.611 | -1.296 to -0.683 |
| ZNF28 | -1.308 | -1.584 to -1.148 | -1.463 to -1.172 |

All six uncentered contrasts are negative. HIVEP2 and GCN1 consequently have different signs from the published signed interaction. That comparison alone is not a contradiction: the publication's reference is the non-targeting guide response, rather than zero change after total-read normalization.

## Disclosed exploratory reference-control diagnostic

After the primary result and the above sign discrepancy were known, the coordinating review requested a non-targeting-reference diagnostic. Its timestamp, outcome knowledge and rules are recorded as a post-primary amendment in `protocol/raw_sensitivity.json`. It does not replace any prespecified estimate. The strict annotation identifies exactly 941 `NO_SITE_` barcodes, matching the publication's stated non-targeting count; exact historical membership equivalence is not independently established.

Their mean hPi-minus-control log2(CPM+1) difference is -0.544245. Subtracting this reference from each gene gives a descriptive difference of differences:

| Gene | NT-centered contrast | Animal omission range | Guide omission range | Published direction agreement |
|---|---:|---:|---:|---|
| HIVEP2 | 0.441 | 0.093 to 0.643 | 0.199 to 0.746 | Yes |
| RNLS | -0.104 | -0.365 to 0.061 | -0.353 to 0.209 | Yes |
| GCN1 | 0.333 | 0.127 to 0.644 | 0.207 to 0.511 | Yes |
| MED11 | -0.248 | -0.462 to 0.172 | -0.486 to 0.001 | Yes |
| ZBED3 | -0.335 | -0.456 to -0.118 | -0.751 to -0.138 | Yes |
| ZNF28 | -0.764 | -1.020 to -0.576 | -0.918 to -0.628 | Yes |

All six centered directions agree with the publication. HIVEP2 and GCN1 remain positive under each omission in this descriptive diagnostic; RNLS and MED11 can cross zero. This explains why reference choice matters. It supplies neither fitted-model p-values nor biological replication, and does not turn nonsignificant published positive scores into supported hits.

## Mapping damage and identity accounting

GCN1 guides map directly to current human GCN1, Entrez 10985, in both original-target and strict chips; no capitalization-based alias guess is used. The mouse library's historical Gcn1l1 label is reconciled to Gcn1 using MGI's unique synonym record.

One intended RNLS guide, `TCCCACACAGCAAGGTACAA`, has current strict matches to both RNLS (55328) and LIPJ (142910). It remains an RNLS *intended* guide in the frozen analysis, while its dual match is flagged. Omitting that guide gives an uncentered contrast of -0.550689 and an NT-centered contrast of approximately -0.006444. Thus it materially affects the magnitude, without changing the negative full-data sign. No unambiguous RNLS-specific mechanistic conclusion follows from that pooled guide alone. The complete annotation audit includes all matches rather than silently choosing one.

## Abundance, guide coverage and mouse-library eligibility

The prespecified matched-background check uses the published signed scores and restricts peers to identical intended guide count plus identical quintiles of mean no-PBMC and pre-transplant CPM. Quintiles are calculated in the complete scored-and-mapped universe. There are 11,474 mapped scored genes; the other 215 retain explicit rows in `results/raw_scored_genes_without_covariates.csv`. All six candidates have available covariates and sufficiently large peer strata.

| Gene | Matched peer count | Probability candidate score exceeds matched peer, with half weight for ties |
|---|---:|---:|
| HIVEP2 | 542 | 0.648 |
| RNLS | 1,695 | 0.337 |
| GCN1 | 346 | 0.942 |
| MED11 | 556 | 0.097 |
| ZBED3 | 433 | 0.062 |
| ZNF28 | 1,214 | 0.127 |

Their unweighted mean is 0.369. This is a descriptive covariate comparison, not the primary theta, a randomization p-value or an adjusted causal effect. It does not suggest that the lack of positive set-level shift is explained by these measured abundance/guide-count covariates. Control abundance is itself post-transplant and may reflect biological selection, so matching cannot remove every selection mechanism.

The mouse screen used GeCKOv2 A. The Addgene `mouse_geckov2_library_a_2.csv` is mapped through the frozen MGI human–mouse homology report, with exact current mouse symbols first and unique MGI historical synonyms second. Human counterparts are never inferred by capitalization. Restricting the published universe to genes with at least one mapped mouse-library member retains all six candidates and 11,104 scored genes; theta is **0.4209467171**, essentially unchanged from the primary 0.4209677880. Homology classes can contain multiple members. This check establishes library availability, not equal guide efficacy, expression or survival eligibility; the absent complete autoimmune count matrix still limits exchangeability assessment.

## Reproduction and output dictionary

Run `python src/fetch_raw_counts.py` to verify/fetch frozen inputs, then `python src/raw_counts.py` after the primary score extraction. The latter uses NumPy and pandas only. `python -m unittest tests.test_raw_counts -v` tests normalization, zero-total rejection, mouse-level contrasts and rank/tie accounting. No fitted count model is used.

* `results/raw_sample_qc.csv`: sample totals and zero-count fractions.
* `results/raw_candidate_guides.csv`: each candidate guide's construct, original gene ID, raw counts and CPM in every graft and baseline sample.
* `results/raw_candidate_summary.csv`: uncentered contrasts, omission ranges and strict-mapping flags.
* `results/raw_leave_one_animal_out.csv` / `raw_leave_one_guide_out.csv`: all prespecified omission estimates.
* `results/raw_exploratory_NT_centered.csv`: explicitly post-primary reference-centered diagnostic and omission ranges.
* `results/raw_matched_background.csv`: same-coverage/abundance peer counts and published-score comparisons.
* `results/raw_mouse_library_sensitivity.json`: restricted-universe coverage and theta.
* `data/processed/raw_annotation_audit.csv`: original intended targets and every strict-sequence match.
* `data/processed/raw_mouse_eligibility.csv`: exact homology class, eligible mouse symbols and historical aliases.

The surviving contribution is a transparent, finite comparison that does not support the predicted positive set shift. Raw checks make its reference choices and biological limitations visible; they do not establish a universal absence of protection, retained glucose responsiveness, a clinical benefit or a new mechanism.

Public annotation authorities are the [Broad GPP CP0043 pool](https://portals.broadinstitute.org/gpp/public/pool/details?pool=CP0043), [Addgene Mouse GeCKOv2 library](https://www.addgene.org/pooled-library/zhang-mouse-gecko-v2/), and [MGI report repository](https://www.informatics.jax.org/downloads/reports/). Cite Doench et al. (2016) for Brunello and Sanjana et al. (2014) for GeCKOv2. Downloading these data does not confer a new redistribution license; source-specific access notes are preserved.
