"""Portable serialization must not weaken the scientific input commitments."""
import csv
from decimal import Decimal
import gzip
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from src.portable_tables import TableSpec, validate_csv_bytes, verify_revision_inputs


ROOT = Path(__file__).resolve().parents[1]
SPEC = TableSpec(
    (("gene", "text"), ("score", "float"), ("q", "float"),
     ("count", "exact"), ("state", "text")),
    "gene", ("score",), ("score",), (("q", .05),))
ROWS = [
    ["gene", "score", "q", "count", "state"],
    ["NA", "0", "0.05", "4", "zero"],
    ["B", "1", "0.01", "4", "positive"],
    ["C", "1", "0.02", "4", "positive"],
    ["D", "3", "0.04", "4", "positive"],
    ["E", "4", "0.06", "4", "positive"],
]


def encode(rows, newline="\r\n"):
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator=newline).writerows(rows)
    return stream.getvalue().encode("utf-8")


class PortableTableTests(unittest.TestCase):
    def test_newlines_numeric_notation_and_sub_femtounit_noise_are_equivalent(self):
        rows = [row.copy() for row in ROWS]
        rows[4][1] = "3.0000000000000004"
        rows[2][2] = "1.0000000000000001e-2"
        rows[2][3] = "4.0"
        changed = encode(rows, "\n")
        self.assertNotEqual(changed, encode(ROWS))
        validate_csv_bytes(changed, encode(ROWS), SPEC)

    def test_material_numeric_and_discrete_changes_are_rejected(self):
        for column, value in [(1, "3.001"), (2, ".04001"), (3, "4.000000000000001")]:
            with self.subTest(column=column):
                rows = [row.copy() for row in ROWS]
                rows[4][column] = value
                with self.assertRaisesRegex(ValueError, "numeric value changed"):
                    validate_csv_bytes(encode(rows), encode(ROWS), SPEC)

    def test_identifiers_text_order_and_completeness_are_exact(self):
        for change in ("renamed", "duplicate", "reorder", "missing", "extra", "label"):
            with self.subTest(change=change):
                rows = [row.copy() for row in ROWS]
                if change == "renamed": rows[1][0] = "NA_changed"
                elif change == "duplicate": rows[1][0] = "B"
                elif change == "reorder": rows[2], rows[3] = rows[3], rows[2]
                elif change == "missing": rows.pop()
                elif change == "extra": rows.append(["F", "4", ".06", "4", "positive"])
                else: rows[2][-1] = "negative"
                with self.assertRaises(ValueError):
                    validate_csv_bytes(encode(rows), encode(ROWS), SPEC)

    def test_schema_changes_missing_values_and_nonfinite_numbers_are_rejected(self):
        for change in ("rename", "reorder", "extra", "width", "", "NaN", "inf", "1e999", "words"):
            with self.subTest(change=change):
                rows = [row.copy() for row in ROWS]
                if change == "rename": rows[0][1] = "different_score"
                elif change == "reorder": rows[0][1], rows[0][2] = rows[0][2], rows[0][1]
                elif change == "extra":
                    for row in rows: row.append("extra")
                elif change == "width": rows[2].pop()
                else: rows[2][1] = change
                with self.assertRaises(ValueError):
                    validate_csv_bytes(encode(rows), encode(ROWS), SPEC)

    def test_tolerance_cannot_change_ties_signs_or_discovery_classification(self):
        for row, column, value, message in [
            (2, 1, "1.0000000000000009", "ranks/ties"),
            (1, 1, "0.0000000000000001", "sign changed"),
            (1, 2, "0.0500000000000001", "decision boundary"),
        ]:
            with self.subTest(message=message):
                rows = [entry.copy() for entry in ROWS]
                rows[row][column] = value
                with self.assertRaisesRegex(ValueError, message):
                    validate_csv_bytes(encode(rows), encode(ROWS), SPEC)

    def test_tolerance_cannot_change_matching_strata(self):
        spec = TableSpec((("gene", "text"), ("covariate", "float"), ("bin", "exact")),
                         "gene", quintile_columns=(("covariate", "bin"),))
        # A tied pair straddles a quintile boundary; its tiny separation changes
        # matching eligibility even though the floating difference is harmless
        # in magnitude. Such a change must fail closed.
        rows = [["gene", "covariate", "bin"]]
        values = [0, 1, 1, 3, 4, 5, 6, 7, 8, 9]
        import pandas as pd
        bins = pd.qcut(values, 5, labels=False, duplicates="drop")
        rows += [[str(i), str(value), str(bin_)] for i, (value, bin_) in enumerate(zip(values, bins))]
        changed = [row.copy() for row in rows]
        changed[3][1] = "1.0000000000000009"
        with self.assertRaisesRegex(ValueError, "matching strata"):
            validate_csv_bytes(encode(changed), encode(rows), spec)

    def test_protocol_and_original_source_hashes_remain_byte_exact(self):
        for filename in ("protocol/analysis_plan.json", "data/raw/original.csv"):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                path = root / filename
                path.parent.mkdir(parents=True)
                source = b'{"unchanged_value": 1}\r\n'
                commitment = {filename: hashlib.sha256(source).hexdigest()}
                path.write_bytes(source)
                verify_revision_inputs(root, commitment)
                path.write_bytes(source.replace(b"\r\n", b"\n"))
                with self.assertRaisesRegex(ValueError, "Revision input changed"):
                    verify_revision_inputs(root, commitment)

    def test_real_regenerated_inputs_are_portable_and_references_cannot_drift(self):
        plan = json.loads((ROOT / "protocol/revision_analysis_plan.json").read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for filename in plan["input_sha256"]:
                target = root / filename
                target.parent.mkdir(parents=True, exist_ok=True)
                data = (ROOT / filename).read_bytes()
                if filename.endswith(".csv"):
                    refname = "protocol/reference_tables/" + Path(filename).name + ".gz"
                    reference = root / refname
                    reference.parent.mkdir(parents=True, exist_ok=True)
                    reference.write_bytes((ROOT / refname).read_bytes())
                    rows = list(csv.reader(io.StringIO(data.decode("utf-8"))))
                    if "reported_q" in rows[0]:
                        col = rows[0].index("reported_q")
                        rows[1][col] = str(Decimal(rows[1][col]) * Decimal("1.0000000000000001"))
                    data = encode(rows, "\n")
                target.write_bytes(data)
            before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
            verify_revision_inputs(root, plan["input_sha256"])
            self.assertTrue(all(p.read_bytes() == data for p, data in before.items()))
            reference = root / "protocol/reference_tables/allogeneic_scores.csv.gz"
            content = gzip.decompress(reference.read_bytes())
            reference.write_bytes(gzip.compress(content.replace(b"CXCL10", b"ALTERED", 1), mtime=0))
            with self.assertRaisesRegex(ValueError, "Frozen table reference changed"):
                verify_revision_inputs(root, plan["input_sha256"])


if __name__ == "__main__":
    unittest.main()
