"""Fetch public raw-count inputs without credentials; freeze source bytes and hashes."""
from pathlib import Path
import concurrent.futures
import datetime
import hashlib
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6008nnn/GSM6008413/suppl/'
FILES = {
    'GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz': BASE + 'GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz',
    'GSM6008413_sample_info.txt.gz': BASE + 'GSM6008413_sample_info.txt.gz',
    'GSM6008413_CP0043_reference_20160120.csv.gz': BASE + 'GSM6008413_CP0043_reference_20160120.csv.gz',
    'GSM6008413.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6008413&targ=self&form=text&view=full',
    'GSE200085.txt': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200085&targ=self&form=text&view=full',
    'brunello_library_contents.txt':'https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt',
    'CP0043_pool.html':'https://portals.broadinstitute.org/gpp/public/pool/details?pool=CP0043',
    'mouse_gecko_addgene.html':'https://www.addgene.org/pooled-library/zhang-mouse-gecko-v2/',
    'CP0043_origtarget_20220324.chip':'https://portals.broadinstitute.org/gpp/public/chipfile/download?filename=CP0043_origtarget_20220324.chip',
    'CP0043_GRCh38_NCBI_CRISPRko_strict_gene_20240130.chip':'https://portals.broadinstitute.org/gpp/public/chipfile/download?filename=CP0043_GRCh38_NCBI_CRISPRko_strict_gene_20240130.chip',
    'mouse_geckov2_library_a_2.csv':'https://media.addgene.org/cms/filer_public/c3/2c/c32ccf3a-00c5-44bd-bca6-eb7aec7ef433/mouse_geckov2_library_a_2.csv',
    'HOM_MouseHumanSequence.rpt':'https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt',
    'MRK_List2.rpt':'https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt',
}

def fetch(item):
    name, url = item
    target = ROOT / 'data/raw/counts' / name
    relative=str(target.relative_to(ROOT)).replace('\\','/')
    manifest=ROOT/'evidence/raw_count_sources.json'
    frozen={r['filename']:r for r in json.loads(manifest.read_text(encoding='utf-8'))} if manifest.exists() else {}
    previous=frozen.get(relative)
    if previous and previous['url']!=url:
        raise ValueError('Frozen source URL changed: '+relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        req = urllib.request.Request(url, headers={'User-Agent':'Andrew3-public-research/1.0'})
        with urllib.request.urlopen(req, timeout=120) as response:
            data=response.read()
        if previous and hashlib.sha256(data).hexdigest()!=previous['sha256']:
            raise ValueError('Remote source drift; no changed input installed: '+relative)
        target.write_bytes(data)
    b = target.read_bytes()
    if previous:
        if hashlib.sha256(b).hexdigest()!=previous['sha256']:
            raise ValueError('Source changed: '+relative)
        return previous
    return {'filename': relative, 'url':url,
            'sha256':hashlib.sha256(b).hexdigest(), 'bytes':len(b),
            'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'access_note':'Open GEO supplementary data / accession record; source attribution retained',
            'role':'Raw-count robustness, sample identity, and barcode-to-construct mapping'}

def main():
    manifest=ROOT/'evidence/raw_count_sources.json'
    old={r['filename']:r for r in json.loads(manifest.read_text(encoding='utf-8'))} if manifest.exists() else {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        records=list(pool.map(fetch, FILES.items()))
    for i,r in enumerate(records):
        if r['filename'] in old:
            if r['sha256']!=old[r['filename']]['sha256']:
                raise ValueError('Source changed; compare and document before accepting: '+r['filename'])
            records[i]=old[r['filename']]
        elif 'broadinstitute' in r['url']:
            r['access_note']='Broad GPP public research-use data; attribution required; source terms retained in CP0043_pool.html'
        elif 'addgene' in r['url']:
            r['access_note']='Public Addgene depositor library annotations; cite Doench et al. 2016 or Sanjana et al. 2014 and Addgene; no new redistribution license claimed'
        elif 'informatics.jax' in r['url']:
            r['access_note']='Public Mouse Genome Informatics report; frozen current snapshot; cite MGI/Jackson Laboratory'
    (ROOT/'evidence').mkdir(exist_ok=True)
    if manifest.exists():
        if json.loads(manifest.read_text(encoding='utf-8'))!=records:
            raise ValueError('Acquisition would alter frozen manifest; explicit versioning required')
    else:
        manifest.write_text(json.dumps(records,indent=2),encoding='utf-8')
    for r in records: print(r['filename'], r['bytes'], r['sha256'])

if __name__=='__main__': main()
