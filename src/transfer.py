"""Frozen, tie-aware finite-table cross-model transfer analysis."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd

def rankdata(values):
    return pd.Series(np.asarray(values)).rank(method="average").to_numpy()

ROOT = Path(__file__).resolve().parents[1]

def superiority(scores, indices):
    scores = np.asarray(scores, dtype=float)
    indices = np.asarray(indices, dtype=int)
    n, m = len(scores), len(indices)
    if not 0 < m < n or len(set(indices)) != m:
        raise ValueError("Need a unique, nonempty proper candidate subset")
    return float((rankdata(scores)[indices].sum() - m * (m + 1) / 2) / (m * (n - m)))

def sample_subsets(rng, n, m, size):
    """Uniform ordered m-tuples, rejecting entire rows with repeated members.

    Conditioning independent uniform tuples on uniqueness produces uniform
    subsets. Row rejection, rather than replacing individual collisions,
    keeps the distribution especially easy to inspect.
    """
    rows = rng.integers(n, size=(size, m))
    while True:
        bad = np.any(np.diff(np.sort(rows, axis=1), axis=1) == 0, axis=1)
        if not np.any(bad):
            return rows
        rows[bad] = rng.integers(n, size=(int(bad.sum()), m))

def competitive_test(scores, indices, draws=1_000_000, seed=20261025):
    ranks = rankdata(scores)
    n, m = len(ranks), len(indices)
    observed_sum = ranks[indices].sum()
    rng = np.random.Generator(np.random.PCG64(seed))
    exceed = 0
    for start in range(0, draws, 50_000):
        rows = sample_subsets(rng, n, m, min(50_000, draws - start))
        exceed += int((ranks[rows].sum(axis=1) >= observed_sum).sum())
    p = (exceed + 1) / (draws + 1)
    return dict(n_universe=n, n_candidates=m, theta=superiority(scores, indices),
                theta_minus_half=superiority(scores, indices)-.5,
                draws=draws, seed=seed, rng="PCG64", null_exceedances=exceed,
                p_one_sided=p, monte_carlo_se=math.sqrt(p*(1-p)/draws))

def load_scores():
    plan = json.loads((ROOT / "protocol/analysis_plan.json").read_text())
    path = ROOT / "data/raw/scores/mmc2.xlsx"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == plan["source_workbook_sha256"]
    d = pd.read_excel(path, sheet_name="Fig3 - in vivo CRISPR screen", header=1)
    assert list(d.columns) == ["gene", "Fold change", "-log10(FDR adj)"]
    assert len(d) == plan["expected_universe_size"] and d.gene.is_unique
    assert d.notna().all().all()
    assert np.isfinite(d.iloc[:, 1:].to_numpy()).all()
    assert (d["-log10(FDR adj)"] >= 0).all()
    d["published_score"] = d["Fold change"]
    d["published_neg_log10_fdr"] = d["-log10(FDR adj)"]
    d["reported_q"] = np.power(10., -d["published_neg_log10_fdr"])
    d["protective_rank"] = rankdata(-d.published_score)
    d["score_percentile"] = (rankdata(d.published_score) - .5) / len(d) * 100
    d["signal_state"] = np.where(d.published_score < 0, "negative",
        np.where(d.published_score == 0, "zero", np.where(d.reported_q <= .05,
        "positive_supported", "positive_uncertain")))
    return d

def map_candidates(d, plan):
    h = pd.read_csv(ROOT / "data/raw/mapping/HOM_MouseHumanSequence.rpt", sep="\t")
    records = []
    for label in plan["original_candidate_labels"]:
        symbol = "Gcn1" if label == "Gcn1l1" else label
        mapped = plan["candidate_mapping"][label]
        item = dict(discovery_label=label, current_mouse_symbol=symbol,
                    human_symbol=mapped, mouse_ncbi_id=None, mouse_mgi_id=None,
                    human_ncbi_id=None, human_hgnc_id=None,
                    mapping_authority="MGI/Alliance HOM_MouseHumanSequence.rpt frozen snapshot",
                    mapping_sensitive=label=="Gm3604", alias_correction="Gcn1l1 -> Gcn1" if label=="Gcn1l1" else "none")
        if mapped:
            mouse = h[(h.Symbol == symbol) & (h["NCBI Taxon ID"] == 10090)]
            assert len(mouse) == 1
            human = h[(h["DB Class Key"] == mouse.iloc[0]["DB Class Key"]) & (h["NCBI Taxon ID"] == 9606)]
            assert len(human) == 1 and human.iloc[0].Symbol == mapped
            item.update(mouse_ncbi_id=int(mouse.iloc[0]["EntrezGene ID"]), mouse_mgi_id=mouse.iloc[0]["Mouse MGI ID"],
                        human_ncbi_id=int(human.iloc[0]["EntrezGene ID"]), human_hgnc_id=human.iloc[0]["HGNC ID"])
        present = mapped in set(d.gene)
        item.update(score_available=present, included=present,
                    exclusion_reason="" if present else ("no_committed_admissible_human_mapping" if not mapped else "absent_reported_score_universe"))
        records.append(item)
    return pd.DataFrame(records)

def main():
    plan = json.loads((ROOT / "protocol/analysis_plan.json").read_text())
    d = load_scores()
    mapping = map_candidates(d, plan)
    genes = plan["primary_human_candidates"]
    assert set(mapping.loc[mapping.included, "human_symbol"]) == set(genes)
    cases = [("primary", genes, d), ("exclude_ZNF28", [g for g in genes if g != "ZNF28"], d)]
    cases += [("leave_out_" + g, [h for h in genes if h != g], d) for g in genes]
    reduced = [g for g in genes if g not in ["RNLS", "HIVEP2"]]
    cases += [("exclude_RNLS_HIVEP2_from_candidates", reduced, d),
              ("exclude_RNLS_HIVEP2_from_universe", reduced, d[~d.gene.isin(["RNLS", "HIVEP2"])].reset_index(drop=True))]
    results = []
    for label, candidates, universe in cases:
        indices = np.flatnonzero(universe.gene.isin(candidates))
        out = competitive_test(universe.published_score.to_numpy(), indices, plan["primary_draws"], plan["seed"])
        out.update(analysis=label, candidate_genes=";".join(candidates),
                   positive_scores=int((universe.iloc[indices].published_score > 0).sum()),
                   positive_supported=int((universe.iloc[indices].signal_state == "positive_supported").sum()),
                   primary_or_sensitivity="primary" if label == "primary" else "sensitivity")
        results.append(out)
    summary = pd.DataFrame(results)
    candidates = d.set_index("gene").loc[genes].reset_index()
    for folder in ("results", "data/processed"):
        (ROOT / folder).mkdir(exist_ok=True, parents=True)
    d.to_csv(ROOT / "data/processed/allogeneic_scores.csv", index=False, float_format="%.17g")
    mapping.to_csv(ROOT / "candidate_mapping.tsv", sep="\t", index=False)
    summary.to_csv(ROOT / "results/transfer_summary.tsv", sep="\t", index=False, float_format="%.17g")
    candidates.to_csv(ROOT / "results/candidate_scores.tsv", sep="\t", index=False, float_format="%.17g")
    pd.concat([summary.assign(record_type="set_test"), candidates.assign(record_type="candidate")], ignore_index=True).to_csv(ROOT / "transfer_results.tsv", sep="\t", index=False, float_format="%.17g")
    audit = dict(workbook_rows=len(d), unique_genes=int(d.gene.nunique()),
        source_sha256=plan["source_workbook_sha256"], discovery_count=11, eligible_count=len(genes),
        coverage_fraction=len(genes)/11, q_encoding="q = 10^(-stored -log10(FDR adj))",
        neg_log10_fdr_distinct=int(d.published_neg_log10_fdr.nunique()),
        neg_log10_fdr_min=float(d.published_neg_log10_fdr.min()),
        neg_log10_fdr_max=float(d.published_neg_log10_fdr.max()),
        noninteger_neg_log10_fdr_values=int((d.published_neg_log10_fdr%1 != 0).sum()),
        tied_score_count=int(d.published_score.duplicated(keep=False).sum()),
        null_sampler="Independent uniform integer tuples conditioned on all entries being distinct; whole duplicate rows resampled",
        score_magnitude_interpretation="Preserve signed published Fold change field; ranks and sign used; no back-transformed biological fold effects",
        inference_scope="competitive gene-label exchangeability only, not biological replicate uncertainty")
    (ROOT / "evidence/score_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(summary.to_string(index=False))
    print(candidates.to_string(index=False))

if __name__ == "__main__":
    main()
