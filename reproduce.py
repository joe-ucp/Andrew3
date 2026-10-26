"""Offline, portable reconstruction from hash-verified public inputs.

Usage: python reproduce.py
Acquire required provider files first: python src/fetch_public_inputs.py
Optional complete historical-source audit: python reproduce.py --source-scope full
Optional manuscript files: python reproduce.py --manuscript
The optional --extra-site-packages directory is APPENDED; it never replaces
already available numerical dependencies. A normal venv needs no such option.
"""
from pathlib import Path
import argparse
import contextlib
import hashlib
import importlib.metadata as metadata
import io
import json
import os
import platform
import runpy
import shutil
import sys
import unittest

ROOT=Path(__file__).resolve().parent

def package_pins():
    return dict(line.strip().split('==',1) for line in (ROOT/'requirements-lock.txt').read_text().splitlines()
                if line.strip() and not line.lstrip().startswith('#'))

def check_environment():
    packages={};errors=[]
    for name,expected in package_pins().items():
        try: actual=metadata.version(name)
        except metadata.PackageNotFoundError: actual='missing'
        packages[name]=actual
        if actual!=expected: errors.append(f'{name}: expected {expected}, found {actual}')
    if errors:
        raise RuntimeError('Pinned environment required. Install requirements-lock.txt in a clean venv.\n'+'\n'.join(errors))
    record={'python_version':platform.python_version(),'python_implementation':platform.python_implementation(),
            'platform':platform.platform(),'packages':packages,
            'numerical_rng':'NumPy Generator(PCG64); primary seed20261025; secondary seeds declared in protocol/revision_analysis_plan.json',
            'network_required_for_reproduction':False,
            'source_scope':os.environ.get('ANDREW3_SOURCE_SCOPE','computational'),
            'environment_note':'Package versions checked; external package directory, if supplied, appended after existing sys.path entries.'}
    (ROOT/'evidence/environment.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

def immutable_inputs():
    # Source bytes, acquisition manifests, curation code and frozen protocols.
    from src.source_inventory import MANIFESTS,inventory
    files=[ROOT/r['filename'] for r in inventory(ROOT)]
    files += [ROOT/p for p in MANIFESTS]
    files += list((ROOT/'protocol').rglob('*'))+[ROOT/'evidence/functional_build.py']
    return {str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files if p.is_file()}

def run_stage(path):
    stream=io.StringIO()
    previous_argv=sys.argv
    try:
        sys.argv=[str(ROOT/path)]
        with contextlib.redirect_stdout(stream): runpy.run_path(str(ROOT/path),run_name='__main__')
    finally:
        sys.argv=previous_argv
    print('Completed '+path)

def hydrate_data_cache(cache):
    """Copy explicit cached raw inputs into a clean checkout; never overwrite drift."""
    cache=cache.resolve();target=(ROOT/'data/raw').resolve()
    if cache==target: return
    if not cache.is_dir(): raise FileNotFoundError('Data cache directory does not exist')
    if not (cache/'scores/mmc2.xlsx').is_file() or not (cache/'literature/manifest.json').is_file():
        raise ValueError('--data-cache must name the retained data/raw directory')
    copied=0
    for source in sorted(cache.rglob('*')):
        if not source.is_file(): continue
        if not source.resolve().is_relative_to(cache):
            raise ValueError('Data cache contains an external file link')
        relative=source.relative_to(cache);destination=(target/relative).resolve()
        if not destination.is_relative_to(target): raise ValueError('Unsafe cached input path')
        if destination.exists():
            if hashlib.sha256(source.read_bytes()).digest()!=hashlib.sha256(destination.read_bytes()).digest():
                raise ValueError('Existing input differs from declared cache: '+str(relative))
            continue
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,destination);copied+=1
    print(f'Restored {copied} raw input files from the declared cache; cache unchanged.')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extra-site-packages',type=Path,help='Optional additional packages directory, appended after existing paths')
    parser.add_argument('--manuscript',action='store_true',help='Also run src/build_manuscript.py; PDF compilation is separately documented')
    parser.add_argument('--data-cache',type=Path,help='Optional retained data/raw directory to copy into a clean export before verifying frozen hashes')
    parser.add_argument('--source-scope',choices=['computational','full'],default='computational',
                        help='computational: verify all original inputs read by the analyses; full: also verify every historical literature/search snapshot')
    args=parser.parse_args()
    if args.extra_site_packages:
        extra=args.extra_site_packages.resolve()
        if not extra.is_dir(): parser.error('Extra package directory does not exist')
        sys.path.append(str(extra))
    os.chdir(ROOT)
    os.environ['ANDREW3_SOURCE_SCOPE']=args.source_scope
    if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
    check_environment()
    if args.data_cache: hydrate_data_cache(args.data_cache)
    from src.source_inventory import write_inventory
    before=immutable_inputs()
    records=write_inventory(ROOT)
    print(f'Verified {len(records)} frozen source files; source scope={args.source_scope}; offline reconstruction.')
    if args.source_scope=='computational':
        print('All analysis stages run. Historical literature/search snapshots are catalogued, not reacquired or re-reviewed by this command.')
    for path in ['src/transfer.py','src/raw_counts.py','src/revision_analysis.py','evidence/functional_build.py','src/figures.py']:
        run_stage(path)
    if args.manuscript:
        builder=ROOT/'src/build_manuscript.py'
        if not builder.exists(): raise FileNotFoundError('--manuscript requested but builder is absent')
        run_stage('src/build_manuscript.py')
    suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'),pattern='test_*.py')
    result=unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful(): raise SystemExit(1)
    if immutable_inputs()!=before:
        raise RuntimeError('Reconstruction mutated a frozen input or manifest')
    print(f'Reconstruction complete: {result.testsRun} tests passed; frozen inputs unchanged.')

if __name__=='__main__': main()
