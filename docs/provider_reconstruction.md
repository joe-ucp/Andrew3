# Reconstruction from a public checkout and provider downloads

The numerical analysis can be reconstructed without the author's retained raw-data cache. The default computational route requires 10 input paths containing nine unique source objects. Every object is checked against the original SHA-256 commitment before use. The complete historical literature/search archive is a separate audit scope; it is not required to rerun the numerical calculations or regenerate the shipped manually curated evidence tables.

## Run the computational route

Use a Python 3.12 environment with the exact packages in `requirements-lock.txt`. From the checkout root, after activating that environment:

```text
python -m pip install -r requirements-lock.txt
python src/fetch_public_inputs.py --fresh
python reproduce.py --source-scope computational
```

`--fresh` proves that none of the computational raw inputs existed at the start of acquisition. For a partially completed download, omit `--fresh` when retrying: every already present input is still checked. The downloader never searches another local directory for a retained cache. It does not alter the source manifests. It exits with code 2 if any required source is blocked, missing or changed, and records individual outcomes in `evidence/provider_acquisition.json`.

The computational reproduction runs every scientific stage: published-score analysis, raw-guide diagnostics, secondary calibration, functional-table serialization, figure generation and the complete test suite. It checks the 10 computational raw paths before and after execution. The generated inventory is `evidence/computational_sources.tsv`; the historical `sources.tsv` is preserved. The console and `evidence/environment.json` state the selected scope.

Add `--manuscript` to regenerate the foundation paper's editable manuscript sources. PDF compilation and the separate `functional_protection/` Perspective have their own documented commands. Numerical reconstruction does not require a TeX installation.

## Downloads and exact input identity

The original hashes and provider URLs are in `evidence/primary_sources.json` and `evidence/raw_count_sources.json`. The downloader reads these commitments directly. The following list is the complete computational source set, not a selection of convenient files.

| Provider object | Installed destination under `data/raw/` | Computational use |
|---|---|---|
| [Sintov supplementary package](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9481918/supplementaryFiles), member `mmc2.xlsx` | `scores/mmc2.xlsx` | Published score universe |
| [MGI homology report](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt) | `mapping/HOM_MouseHumanSequence.rpt` and `counts/HOM_MouseHumanSequence.rpt` | Candidate mapping and mouse-library eligibility |
| [GEO count matrix](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6008nnn/GSM6008413/suppl/GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz) | `counts/GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz` | Raw guide counts |
| [GEO sample information](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6008nnn/GSM6008413/suppl/GSM6008413_sample_info.txt.gz) | `counts/GSM6008413_sample_info.txt.gz` | Sample allocation |
| [GEO barcode reference](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6008nnn/GSM6008413/suppl/GSM6008413_CP0043_reference_20160120.csv.gz) | `counts/GSM6008413_CP0043_reference_20160120.csv.gz` | Barcode/construct mapping |
| [Broad original targets](https://portals.broadinstitute.org/gpp/public/chipfile/download?filename=CP0043_origtarget_20220324.chip) | `counts/CP0043_origtarget_20220324.chip` | Original guide targets |
| [Broad strict annotation](https://portals.broadinstitute.org/gpp/public/chipfile/download?filename=CP0043_GRCh38_NCBI_CRISPRko_strict_gene_20240130.chip) | `counts/CP0043_GRCh38_NCBI_CRISPRko_strict_gene_20240130.chip` | Guide identity sensitivity |
| [Addgene mouse library](https://media.addgene.org/cms/filer_public/c3/2c/c32ccf3a-00c5-44bd-bca6-eb7aec7ef433/mouse_geckov2_library_a_2.csv) | `counts/mouse_geckov2_library_a_2.csv` | Discovery-library eligibility |
| [MGI marker report](https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt) | `counts/MRK_List2.rpt` | Historical mouse aliases |

The supplementary ZIP is a transport container. The computational commitment is the exact extracted workbook. Providers can rebuild a ZIP without changing that workbook. The downloader records the downloaded container's hash and verifies the workbook's original hash; it does not replace the historical ZIP commitment or claim to have recovered that historical container.

## When a provider blocks automated download

Use the listed public provider links in a browser or the provider's documented download service. Save the original binary files, retaining the names above, into a separate download folder. Keep the GEO files compressed as `.gz`. For the supplementary package, extract `mmc2.xlsx`. Then run:

```text
python src/fetch_public_inputs.py --manual-downloads PATH_TO_PROVIDER_DOWNLOADS
python reproduce.py --source-scope computational
```

The manual importer accepts either the basename or the complete project-relative path inside that folder. It verifies the same original hashes before copying anything and does not modify the downloaded originals. A download failure is not permission to bypass an access restriction. A hash mismatch is not an acceptable replacement: obtain the original dated bytes or report the reconstruction as blocked. In particular, MGI reports are live reports and can change. If a provider stops serving the committed version, the exact historical analysis requires a lawful versioned deposit of those bytes; changing a hash would create a different analysis.

## Historical evidence audit is a different scope

```text
python reproduce.py --source-scope full --data-cache PATH_TO_RETAINED_DATA_RAW
```

This optional route verifies all historical payload commitments in the seven manifests before running the same stages and tests. It requires the retained historical files. The computational provider downloader does not claim to reproduce old database result sets, browser text extractions, rendered figure pages, or manual literature judgments. Those historical searches can change, some source accesses were restricted, and derived screenshots/text are not necessarily served as identical bytes by a provider.

The public checkout supplies the search protocols and logs, comparison-level curation, source identifiers, original hash catalog, and the deterministic table builder. Rebuilding those tables checks serialization and scientific contracts; it is not an independent re-extraction or validation of the literature evidence. The full historical mode remains available and continues to fail on missing or changed committed payloads. No source checks are silently bypassed under the name of full reconstruction.

In computational mode, each curated source-file locator must resolve to both an original manifest commitment and a historical catalog entry with a nonempty size, SHA-256 and source provenance. Unknown locators fail. If a historical payload is present, its bytes are also checked; otherwise its status is explicitly commitment-only, not payload-verified. Full mode requires the physical payload as well.

## Executed acquisition check

On 26 October 2026 UTC, a test started with no computational input files and downloaded all nine unique objects directly from the providers. All 10 installed paths matched their original hashes; no manual download or retained cache was used. Both source manifests were unchanged. The acquisition receipt is `evidence/provider_acquisition_test.json`.

The Europe PMC ZIP container differed from the historical container hash, while its `mmc2.xlsx` member matched exactly. That distinction is recorded rather than treated as source drift in the workbook. Current successful acquisition does not guarantee perpetual provider availability. The final release reconstruction receipt separately records the downstream numerical and test results.
