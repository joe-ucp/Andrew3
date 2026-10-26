"""Acquire the frozen v3 uniform literature searches; never overwrite raw inputs.

Offline reconstruction consumes the retained snapshots and does not run this helper.
Protocol and authoritative metadata are frozen before publication queries.
"""
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import datetime
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/raw/literature/uniform20261025'
PLAN = ROOT / 'protocol/functional_search_v3.json'
EXPECTED = '951e7160b4f1f180d2e22f95ed16643a7084c98b8a9a8fb35c93ebfcdaf97268'
sys.stdout.reconfigure(encoding='utf-8')
assert hashlib.sha256(PLAN.read_bytes()).hexdigest() == EXPECTED
plan = json.loads(PLAN.read_text(encoding='utf-8'))
if list(OUT.glob('*_uniform_search.json')):
    raise RuntimeError('Uniform search snapshot exists; do not overwrite or mix acquisitions.')
manifest = json.loads((OUT / 'alias_manifest.json').read_text())
for entry in plan['targets']:
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urlencode(
        dict(query=entry['exact_query'], format='json', resultType='core', pageSize=1000))
    path = OUT / (entry['target'] + '_uniform_search.json')
    assert not path.exists()
    with urlopen(Request(url, headers={'User-Agent':'literature bounded evidence assessment'}), timeout=90) as response:
        data = response.read()
        final_url = response.url
    decoded = json.loads(data)
    assert 'hitCount' in decoded and 'resultList' in decoded
    path.write_bytes(data)
    manifest.append(dict(filename=path.relative_to(ROOT).as_posix(), url=url,
        final_url=final_url, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
        retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        role='Frozen uniform alias-aware title/abstract search; target ' + entry['target'],
        access_note='Public Europe PMC metadata; original source rights retained'))
    (OUT / 'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    records = decoded.get('resultList',{}).get('result',[])
    print(entry['target'], 'hits', decoded['hitCount'], 'returned', len(records), flush=True)
    time.sleep(.2)
