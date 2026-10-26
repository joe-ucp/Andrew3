"""Validate regenerated tables without treating CSV serialization as science.

Original inputs and protocols retain exact SHA-256 commitments. Only the three
declared, regenerated CSVs permit different newlines/numeric representations.
Their reference CSV bytes are compressed, unchanged copies of the original
commitments, never learned or updated during reproduction.
"""
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
import csv
import gzip
import hashlib
import io
import math
import re

import numpy as np
import pandas as pd


# Reference-relative floating error only; discrete values and scientific
# ordering/decision boundaries are checked separately and remain exact.
FLOAT_RTOL = Decimal("1e-13")
FLOAT_ATOL = Decimal("1e-15")
NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")


@dataclass(frozen=True)
class TableSpec:
    columns: tuple
    identifier: str
    rank_columns: tuple = ()
    sign_columns: tuple = ()
    cutoff_columns: tuple = ()
    quintile_columns: tuple = ()


SCORE_COLUMNS = (
    ("gene", "text"), ("Fold change", "float"), ("-log10(FDR adj)", "float"),
    ("published_score", "float"), ("published_neg_log10_fdr", "float"),
    ("reported_q", "float"), ("protective_rank", "exact"),
    ("score_percentile", "float"), ("signal_state", "text"),
)
TABLES = {
    "data/processed/allogeneic_scores.csv": TableSpec(
        SCORE_COLUMNS, "gene", ("published_score",), ("published_score",),
        (("reported_q", .05),)),
    "data/processed/raw_scored_universe_covariates.csv": TableSpec(
        SCORE_COLUMNS + (("guide_count", "exact"), ("mean_control_CPM", "float"),
        ("mean_hPi_CPM", "float"), ("mean_pretransplant_CPM", "float"),
        ("mean_library_CPM", "float"), ("descriptive_log2_CPM_contrast", "float"),
        ("control_quintile", "exact"), ("pretransplant_quintile", "exact")),
        "gene", ("published_score", "descriptive_log2_CPM_contrast"),
        ("published_score", "descriptive_log2_CPM_contrast"), (("reported_q", .05),),
        (("mean_control_CPM", "control_quintile"),
         ("mean_pretransplant_CPM", "pretransplant_quintile"))),
    "results/raw_exploratory_NT_sample_reference.csv": TableSpec(
        (("sample", "text"), ("NT_mean_log2_CPM_plus1", "float")), "sample"),
}


def _rows(data, spec, label):
    try:
        rows = list(csv.reader(io.StringIO(data.decode("utf-8"), newline=""), strict=True))
    except (UnicodeError, csv.Error) as exc:
        raise ValueError(f"{label}: invalid UTF-8 CSV") from exc
    columns = [column for column, _ in spec.columns]
    if not rows or rows[0] != columns:
        raise ValueError(f"{label}: schema changed (column names/order)")
    rows = rows[1:]
    if any(len(row) != len(columns) for row in rows):
        raise ValueError(f"{label}: malformed row width")
    identifiers = [row[columns.index(spec.identifier)] for row in rows]
    if any(not key for key in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError(f"{label}: missing or duplicate identifier")
    return rows


def _number(value, label):
    if not NUMBER.fullmatch(value):
        raise ValueError(f"{label}: missing, nonfinite or invalid numeric value")
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f"{label}: invalid numeric value") from exc
    if not number.is_finite() or not math.isfinite(float(number)):
        raise ValueError(f"{label}: nonfinite numeric value")
    return number


def validate_csv_bytes(actual, reference, spec, label="regenerated table"):
    """Require identical ordered records/schema and bounded numerical changes.

    Row order matters because it determines seeded sampling. Identifiers and
    text labels are literal strings (including strings such as NA). No sorting,
    coercion of missing values, or column intersection hides changed content.
    """
    observed = _rows(actual, spec, label)
    expected = _rows(reference, spec, label + " reference")
    if len(observed) != len(expected):
        raise ValueError(f"{label}: row count changed")
    for row_number, (row, baseline) in enumerate(zip(observed, expected), start=2):
        for (column, kind), value, old in zip(spec.columns, row, baseline):
            location = f"{label}: row {row_number}, {column}"
            if kind == "text":
                if value != old:
                    raise ValueError(f"{location}: identifier/text/order changed")
                continue
            number, original = _number(value, location), _number(old, location)
            tolerance = (max(FLOAT_ATOL, FLOAT_RTOL * abs(original))
                         if kind == "float" else Decimal(0))
            if abs(number - original) > tolerance:
                raise ValueError(f"{location}: numeric value changed beyond tolerance")

    # Match the downstream pandas parser as well as decimal serialization.
    # Tiny perturbations may pass a numerical bound but must not change ranks,
    # ties, sign, discovery classification, or the frozen matching strata.
    current = pd.read_csv(io.BytesIO(actual))
    baseline = pd.read_csv(io.BytesIO(reference))
    for column in spec.rank_columns:
        if not np.array_equal(current[column].rank(), baseline[column].rank()):
            raise ValueError(f"{label}: ranks/ties changed in {column}")
    for column in spec.sign_columns:
        if not np.array_equal(np.sign(current[column]), np.sign(baseline[column])):
            raise ValueError(f"{label}: sign changed in {column}")
    for column, cutoff in spec.cutoff_columns:
        if not np.array_equal(current[column] <= cutoff, baseline[column] <= cutoff):
            raise ValueError(f"{label}: decision boundary changed in {column}")
    for covariate, stratum in spec.quintile_columns:
        bins = pd.qcut(current[covariate], 5, labels=False, duplicates="drop")
        if not np.array_equal(bins, baseline[stratum]):
            raise ValueError(f"{label}: matching strata changed for {covariate}")


def verify_revision_inputs(root, commitments):
    """Validate every frozen commitment; never rewrite an input or reference."""
    root = Path(root)
    for filename, expected_hash in commitments.items():
        actual = (root / filename).read_bytes()
        if filename not in TABLES:
            if hashlib.sha256(actual).hexdigest() != expected_hash:
                raise ValueError("Revision input changed: " + filename)
            continue
        reference_path = root / "protocol/reference_tables" / (Path(filename).name + ".gz")
        reference = gzip.decompress(reference_path.read_bytes())
        if hashlib.sha256(reference).hexdigest() != expected_hash:
            raise ValueError("Frozen table reference changed: " + filename)
        validate_csv_bytes(actual, reference, TABLES[filename], filename)
