"""Acquire computational inputs from providers, retaining original SHA-256 checks.

Run from a public checkout: python src/fetch_public_inputs.py --fresh
Blocked downloads can be supplied with --manual-downloads DIRECTORY. No retained
local cache is consulted, no source commitment is rewritten, and changed bytes
are rejected before installation. See docs/provider_reconstruction.md.
"""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import hashlib
import io
import json
import sys
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from src.source_inventory import COMPUTATIONAL_INPUTS,manifest_records,relative_path

def digest(data):
    return hashlib.sha256(data).hexdigest()

def checked_write(root,relative,data,expected):
    """Never install unverified bytes or overwrite a conflicting local file."""
    relative=relative_path(relative)
    target=(root/relative).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError('Input path escapes checkout: '+relative)
    actual=digest(data)
    if actual!=expected:
        raise ValueError(f'Provider/source drift: expected {expected}, received {actual}')
    if target.exists():
        if digest(target.read_bytes())!=expected:
            raise ValueError('Existing local input differs from frozen commitment: '+relative)
        return
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(data)

def read_provider(url,timeout):
    req=urllib.request.Request(url,headers={'User-Agent':'Andrew3-public-research/2.0'})
    with urllib.request.urlopen(req,timeout=timeout) as response:
        return response.read()

def acquire_one(root,records,manual=None,timeout=120):
    """One byte-identical source may satisfy more than one declared input path."""
    source=records[0]
    relative=source['filename'];expected=source['sha256']
    details={'url':source['url'],'expected_sha256':expected}
    data=None;status=None
    for r in records:
        target=root/r['filename']
        if target.exists():
            data=target.read_bytes()
            if digest(data)!=expected:
                raise ValueError('Existing local input differs from frozen commitment: '+r['filename'])
            status='verified_existing';break
    if data is None and manual:
        candidates=[manual/relative,manual/Path(relative).name]
        for path in candidates:
            if path.is_file():
                data=path.read_bytes();status='imported_manual_verified';break
    if data is None:
        data=read_provider(source['url'],timeout)
        status='downloaded_verified'
        if source.get('archive_member'):
            details['download_container_sha256']=digest(data)
            with zipfile.ZipFile(io.BytesIO(data)) as archive:
                member=source['archive_member']
                if member not in archive.namelist():
                    matches=[n for n in archive.namelist() if n.endswith('/'+member)]
                    if len(matches)!=1:
                        raise ValueError('Expected workbook member is missing or ambiguous')
                    member=matches[0]
                data=archive.read(member)
                details['archive_member']=member
            status='downloaded_archive_member_verified'
    if digest(data)!=expected:
        raise ValueError(f'Provider/source drift: expected {expected}, received {digest(data)}; no changed input installed')
    for r in records:
        checked_write(root,r['filename'],data,expected)
    return [dict(filename=r['filename'],status=status,bytes=len(data),**details) for r in records]

def acquire(root=ROOT,manual=None,timeout=120,fresh=False,workers=3):
    root=root.resolve()
    commitments=manifest_records(root,'computational')
    records=[dict(r,filename=relative_path(r['filename'])) for _,r in commitments]
    before={name:digest((root/name).read_bytes()) for name in {m for m,_ in commitments}}
    if fresh and any((root/r['filename']).exists() for r in records):
        raise ValueError('--fresh requires all computational raw inputs to be absent before acquisition')
    groups={}
    for record in records:
        key=(record['url'],record['sha256'],record.get('archive_member',''))
        groups.setdefault(key,[]).append(record)
    outcomes=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(acquire_one,root,rs,manual,timeout):rs for rs in groups.values()}
        for future in concurrent.futures.as_completed(futures):
            try: outcomes.extend(future.result())
            except Exception as exc:
                outcomes.extend(dict(filename=r['filename'],url=r['url'],expected_sha256=r['sha256'],
                    status='unavailable_or_changed',error=f'{type(exc).__name__}: {exc}') for r in futures[future])
    outcomes.sort(key=lambda r:r['filename'])
    after={name:digest((root/name).read_bytes()) for name in before}
    if before!=after:
        raise RuntimeError('Source acquisition mutated a frozen manifest')
    ok=all(r['status']!='unavailable_or_changed' for r in outcomes)
    report={'schema_version':1,'scope':'computational','status':'PASS' if ok else 'INCOMPLETE',
        'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'fresh_inputs_required':fresh,'retained_cache_used':False,
        'source_paths_required':len(COMPUTATIONAL_INPUTS),'unique_download_objects':len(groups),
        'frozen_manifests_unchanged':True,'manifest_sha256':before,'files':outcomes,
        'historical_evidence_scope':'Historical literature/search snapshots were not acquired or verified. Shipped manual curation is a computational input, not an automated literature re-extraction.',
        'manual_resolution':'Download missing files from the listed provider URLs into a separate folder using the destination basename, then rerun with --manual-downloads. Hash mismatch must be resolved with the original dated bytes; current substitute bytes cannot reconstruct this frozen analysis.'}
    (root/'evidence').mkdir(exist_ok=True)
    (root/'evidence/provider_acquisition.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT,help='Public checkout containing frozen source manifests')
    parser.add_argument('--manual-downloads',type=Path,help='Provider downloads by basename or project-relative path; verified before import')
    parser.add_argument('--fresh',action='store_true',help='Require every computational raw input to be initially absent')
    parser.add_argument('--timeout',type=int,default=120)
    args=parser.parse_args()
    report=acquire(args.root,args.manual_downloads,args.timeout,args.fresh)
    print(json.dumps(report,indent=2))
    if report['status']!='PASS': raise SystemExit(2)

if __name__=='__main__':main()
