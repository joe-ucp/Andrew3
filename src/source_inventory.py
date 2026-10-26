"""Verify retained source commitments and merge them without changing manifests."""
from pathlib import Path, PurePosixPath
import csv
import hashlib
import json
import os
import zipfile

ROOT=Path(__file__).resolve().parents[1]
MANIFESTS=['evidence/primary_sources.json','evidence/raw_count_sources.json',
           'data/raw/literature/manifest.json',
           'data/raw/literature/revision20261025/manifest.json',
           'data/raw/literature/uniform20261025/manifest.json',
           'data/raw/literature/revision3_citations/manifest.json',
           'data/raw/literature/bode2023_final/manifest.json']
FIELDS=['url_or_accession','retrieval_utc','retrieval_basis','filename','sha256','bytes',
        'access_note','role','source_status','source_manifest','parent_source','archive_member']

# These are the raw bytes actually opened by transfer.py and raw_counts.py.
# Later numerical stages consume their outputs. Functional tables are rebuilt
# from the shipped, manually adjudicated evidence/functional_build.py records;
# regenerating those tables does not repeat the historical literature search.
COMPUTATIONAL_INPUTS={
    'data/raw/scores/mmc2.xlsx': 'Published-score analysis',
    'data/raw/mapping/HOM_MouseHumanSequence.rpt': 'Candidate mapping',
    'data/raw/counts/GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz': 'Raw-guide counts',
    'data/raw/counts/GSM6008413_sample_info.txt.gz': 'Sample allocation',
    'data/raw/counts/GSM6008413_CP0043_reference_20160120.csv.gz': 'Barcode/construct reference',
    'data/raw/counts/CP0043_origtarget_20220324.chip': 'Original guide targets',
    'data/raw/counts/CP0043_GRCh38_NCBI_CRISPRko_strict_gene_20240130.chip': 'Strict guide annotation',
    'data/raw/counts/mouse_geckov2_library_a_2.csv': 'Mouse-library eligibility',
    'data/raw/counts/HOM_MouseHumanSequence.rpt': 'Mouse-library homology',
    'data/raw/counts/MRK_List2.rpt': 'Mouse-library historical aliases',
}

def source_scope(scope=None):
    value=scope or os.environ.get('ANDREW3_SOURCE_SCOPE','computational')
    if value not in ('computational','full'):
        raise ValueError('Source scope must be computational or full')
    return value

def manifest_records(root=ROOT,scope=None):
    """Read original commitments without requiring retained payloads."""
    scope=source_scope(scope)
    names=MANIFESTS if scope=='full' else MANIFESTS[:2]
    records=[]
    for name in names:
        manifest=json.loads((root/name).read_text(encoding='utf-8'))
        items=manifest['files'] if isinstance(manifest,dict) else manifest
        for record in items:
            filename=relative_path(record.get('filename',record.get('local')))
            if scope=='full' or filename in COMPUTATIONAL_INPUTS:
                records.append((name,record))
    if scope=='computational':
        found={relative_path(r.get('filename',r.get('local'))) for _,r in records}
        if found!=set(COMPUTATIONAL_INPUTS):
            raise ValueError('Computational source commitments are incomplete: '+str(sorted(set(COMPUTATIONAL_INPUTS)-found)))
    return records

def relative_path(value):
    p=PurePosixPath(str(value).replace('\\','/'))
    if p.is_absolute() or '..' in p.parts or ':' in str(p):
        raise ValueError('Expected a portable project-relative source path: '+str(p))
    return str(p)

def inventory(root=ROOT,scope=None):
    scope=source_scope(scope)
    rows=[]
    for manifest_name in (MANIFESTS if scope=='full' else MANIFESTS[:2]):
        manifest=json.loads((root/manifest_name).read_text(encoding='utf-8'))
        items=manifest['files'] if isinstance(manifest,dict) else manifest
        original={r.get('source_id'):r for r in items if r.get('url') and r.get('retrieved_utc')}
        for r in items:
            name=relative_path(r.get('filename',r.get('local')))
            if scope=='computational' and name not in COMPUTATIONAL_INPUTS:
                continue
            p=root/name
            if not p.is_file(): raise FileNotFoundError('Missing frozen source: '+name)
            b=p.read_bytes();digest=hashlib.sha256(b).hexdigest()
            if digest!=r['sha256']: raise ValueError('Frozen source hash mismatch: '+name)
            parent=r.get('parent_source','')
            origin=original.get(parent,{})
            # PDF/archive-member derivatives inherit acquisition from the package.
            if parent and len(PurePosixPath(name).parts)>4:
                options=[v for k,v in original.items() if k and k.startswith(parent+'_supplementary')]
                if options: origin=options[0]
            url=r.get('url',origin.get('url',''))
            date=r.get('retrieved_utc',origin.get('retrieved_utc',''))
            basis='recorded_source_acquisition'
            if not r.get('retrieved_utc') and origin: basis='inherited_parent_source_acquisition'
            if isinstance(manifest,dict):
                date=manifest.get('retrieval_date_utc','');basis='recorded_date_only'
            if not url or not date: raise ValueError('Unresolved URL/date provenance for '+name)
            role=r.get('role') or r.get('derivation','Primary literature for candidate-specific evidence')
            status='retained'
            if len(b)==0:
                status='empty_asset_not_evidence';role='Empty source asset; unavailable for evidence extraction'
            elif p.suffix.lower()=='.zip':
                with zipfile.ZipFile(p) as z:
                    if not any(v.file_size for v in z.infolist()):
                        status='archive_only_empty_assets';role='Archive retrieval yielded only empty asset; not evidence'
            rows.append(dict(url_or_accession=url,retrieval_utc=date,retrieval_basis=basis,
                filename=name,sha256=digest,bytes=len(b),
                access_note=r.get('access_note',r.get('access_license_note','Public repository snapshot; original rights retained')),
                role=role,source_status=status,source_manifest=manifest_name,
                parent_source=parent,archive_member=r.get('archive_member','')))
    by_name={}
    for r in rows:
        if r['filename'] in by_name and by_name[r['filename']]['sha256']!=r['sha256']:
            raise ValueError('Conflicting source commitments: '+r['filename'])
        by_name[r['filename']]=r
    if scope=='computational' and set(by_name)!=set(COMPUTATIONAL_INPUTS):
        raise ValueError('Computational source commitments are incomplete')
    return sorted(by_name.values(),key=lambda r:r['filename'])

def write_inventory(root=ROOT,scope=None):
    scope=source_scope(scope)
    rows=inventory(root,scope)
    with inventory_path(root,scope).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    return rows

def inventory_path(root=ROOT,scope=None):
    return root/('sources.tsv' if source_scope(scope)=='full' else 'evidence/computational_sources.tsv')

def validate_source_locator(root,filename,scope=None):
    """Check a curated source pointer without claiming unfetched bytes verified.

    Computational mode still requires an exact original manifest commitment,
    nonempty historical size, and resolved provenance in the historical catalog.
    Any present payload must match; full mode additionally requires its presence.
    """
    filename=relative_path(filename);scope=source_scope(scope)
    matches=[(m,r) for m,r in manifest_records(root,'full')
             if relative_path(r.get('filename',r.get('local')))==filename]
    if not matches:
        raise ValueError('Unknown historical source locator: '+filename)
    with (root/'sources.tsv').open(encoding='utf-8',newline='') as stream:
        rows=[r for r in csv.DictReader(stream,delimiter='\t') if r['filename']==filename]
    if len(rows)!=1:
        raise ValueError('Historical source locator has no unique catalog entry: '+filename)
    row=rows[0];expected=row['sha256']
    if len(expected)!=64 or any(c not in '0123456789abcdef' for c in expected) or int(row['bytes'])<=0:
        raise ValueError('Historical source lacks a nonempty size/hash commitment: '+filename)
    for manifest_name,record in matches:
        if record['sha256']!=expected or row['source_manifest']!=manifest_name:
            raise ValueError('Historical catalog/manifest commitment mismatch: '+filename)
        if record.get('bytes') is not None and int(record['bytes'])!=int(row['bytes']):
            raise ValueError('Historical catalog/manifest size mismatch: '+filename)
    if not all(row.get(k) for k in ('url_or_accession','retrieval_utc','retrieval_basis','access_note','role')):
        raise ValueError('Historical source locator lacks resolved provenance: '+filename)
    path=root/filename
    if path.is_file():
        data=path.read_bytes()
        if len(data)!=int(row['bytes']) or hashlib.sha256(data).hexdigest()!=expected:
            raise ValueError('Historical source payload differs from commitment: '+filename)
        status='payload_hash_and_size_verified'
    else:
        if scope=='full':raise FileNotFoundError('Missing historical source payload: '+filename)
        status='historical_commitment_only_payload_not_fetched'
    return dict(filename=filename,status=status,sha256=expected,bytes=int(row['bytes']))

if __name__=='__main__':
    records=write_inventory();print(f'Verified {len(records)} source files; explicit scope: {source_scope()}.')
