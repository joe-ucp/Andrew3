"""Retrieve immutable input bytes; fail on changed primary workbook."""
from pathlib import Path
import datetime
import hashlib
import io
import json
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def fetch(url, relative, role, expected=None):
    path = ROOT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    manifest = ROOT / "evidence/primary_sources.json"
    frozen = {}
    if manifest.exists():
        frozen = {r["filename"]: r for r in json.loads(manifest.read_text(encoding="utf-8"))}
    previous = frozen.get(relative)
    if previous:
        if previous["url"] != url:
            raise ValueError(f"Frozen source URL changed: {relative}")
        if expected and expected != previous["sha256"]:
            raise ValueError(f"Conflicting source commitments: {relative}")
        expected = previous["sha256"]
    if not path.exists():
        req = urllib.request.Request(url, headers={"User-Agent": "Andrew3-public-research/1.0"})
        with urllib.request.urlopen(req, timeout=120) as response:
            data = response.read()
        if expected and hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f"Remote source drift; frozen input not replaced: {relative}")
        path.write_bytes(data)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if expected and digest != expected:
        raise ValueError(f"Source hash changed: {relative}: {digest}")
    if previous:
        return previous
    return dict(url=url, filename=relative, sha256=digest,
                retrieved_utc=datetime.datetime.fromtimestamp(path.stat().st_mtime, datetime.timezone.utc).isoformat(),
                role=role, access_note="Public source; third-party rights retained")

def main():
    records = []
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9481918/supplementaryFiles"
    records.append(fetch(url, "data/raw/scores/PMC9481918_supplementary.zip", "Sintov original supplementary package"))
    with zipfile.ZipFile(ROOT / records[-1]["filename"]) as archive:
        matches = [n for n in archive.namelist() if n.endswith("mmc2.xlsx")]
        if len(matches) != 1:
            raise ValueError(f"Ambiguous workbook members: {archive.namelist()}")
        data = archive.read(matches[0])
    path = ROOT / "data/raw/scores/mmc2.xlsx"
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError("Retained workbook differs from frozen archive member")
    else:
        if hashlib.sha256(data).hexdigest() != "d08660be0403e60f7ea8f2f025ae1d922862d4c4922fb5c0a6e8a3c9ac0920a8":
            raise ValueError("Archive workbook member no longer matches source commitment")
        path.write_bytes(data)
    record = fetch(url, "data/raw/scores/mmc2.xlsx", "Published allogeneic score universe", "d08660be0403e60f7ea8f2f025ae1d922862d4c4922fb5c0a6e8a3c9ac0920a8")
    record["archive_member"] = matches[0]
    records.append(record)
    records.append(fetch("https://www.informatics.jax.org/marker/MGI:3781781", "data/raw/mapping/Gm3604_MGI.html", "Mapping-sensitive Gm3604 to ZNF28 authority"))
    records.append(fetch("https://www.informatics.jax.org/marker/MGI:2444248", "data/raw/mapping/Gcn1_MGI.html", "Gcn1l1 alias and GCN1 orthology authority"))
    records.append(fetch("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", "data/raw/mapping/HOM_MouseHumanSequence.rpt", "MGI mouse-human homology snapshot"))
    manifest = ROOT / "evidence/primary_sources.json"
    if manifest.exists():
        if json.loads(manifest.read_text(encoding="utf-8")) != records:
            raise ValueError("Acquisition would alter frozen source manifest; explicit versioning required")
    else:
        manifest.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(json.dumps([{k: r[k] for k in ("filename", "sha256")} for r in records], indent=2))

if __name__ == "__main__":
    main()
