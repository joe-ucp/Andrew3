"""One-time acquisition of the bounded 2026-10-25 revision search.

Offline reproduction uses the retained cache. This acquisition helper refuses
to replace an existing snapshot or partial acquisition; a new search requires
an explicitly versioned output location and protocol amendment.
"""
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import json, hashlib, datetime, time, sys
sys.stdout.reconfigure(encoding='utf-8')

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/raw/literature/revision20261025'
if OUT.exists() and any(OUT.iterdir()):
    raise RuntimeError('Revision source snapshot already exists; no files changed. Use the retained cache for reproduction. A fresh acquisition requires an explicitly versioned location.')
OUT.mkdir(parents=True,exist_ok=True)
queries={
 'LGSN':'(TITLE_ABS:LGSN OR TITLE_ABS:lengsin OR TITLE_ABS:GLULD1) AND (TITLE_ABS:insulin OR TITLE_ABS:islet OR TITLE_ABS:beta-cell)',
 'PLIN4':'(TITLE_ABS:PLIN4 OR TITLE_ABS:"perilipin 4" OR TITLE_ABS:"S3-12") AND (TITLE_ABS:insulin OR TITLE_ABS:islet OR TITLE_ABS:beta-cell)',
 'TGM6':'(TITLE_ABS:TGM6 OR TITLE_ABS:"transglutaminase 6" OR TITLE_ABS:TGM3L) AND (TITLE_ABS:insulin OR TITLE_ABS:islet OR TITLE_ABS:beta-cell)',
 'GRM2':'(TITLE_ABS:GRM2 OR TITLE_ABS:mGluR2 OR TITLE_ABS:mGlu2 OR TITLE_ABS:GPRC1B) AND (TITLE_ABS:insulin OR TITLE_ABS:islet OR TITLE_ABS:beta-cell)',
 'Or8b47':'(TITLE_ABS:Olfr909 OR TITLE_ABS:Or8b47 OR TITLE_ABS:Olfr911 OR TITLE_ABS:"MOR165-1" OR TITLE_ABS:"MOR166-1") AND (TITLE_ABS:insulin OR TITLE_ABS:islet OR TITLE_ABS:beta-cell)',
 'Brice2002':'EXT_ID:11935156 AND SRC:MED',
 'TGM6context2021':'DOI:10.1038/s41598-021-98792-1',
 'PLIN5context2015':'PMCID:PMC4375085',
 'PLIN4context2026':'EXT_ID:41391763 AND SRC:MED',
 'PLIN4context2013':'PMCID:PMC3625749',
 'PRISMA2020':'DOI:10.1136/bmj.n71',
 'PRISMAS2021':'DOI:10.1186/s13643-020-01542-z',
 'Lakens2017':'DOI:10.1177/1948550617697177',
}
manifest=[]
def save(name,url):
    req=Request(url,headers={'User-Agent':'literature evidence assessment; source snapshot','Accept':'*/*'})
    with urlopen(req,timeout=60) as r: data=r.read(); final=r.url
    path=OUT/name;path.write_bytes(data)
    manifest.append(dict(filename=str(path.relative_to(ROOT)).replace('\\','/'),url=url,final_url=final,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),role='Revision literature search, primary contextual/method source or authoritative gene metadata',access_note='Publicly accessible source snapshot; original rights retained; not new biological validation'))
    return data
for name,q in queries.items():
    q='('+q+') AND FIRST_PDATE:[* TO 2026-10-25]'
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode(dict(query=q,format='json',resultType='core',pageSize=1000))
    data=json.loads(save(name+'_search.json',url))
    hits=data.get('hitCount',0)
    print(name,hits,[(x.get('id'),x.get('title','')[:160]) for x in data.get('resultList',{}).get('result',[])][:15])
    if name in ['TGM6context2021','PLIN4context2026','PLIN4context2013','PLIN5context2015','PRISMA2020','PRISMAS2021','Lakens2017']:
        results=data.get('resultList',{}).get('result',[])
        if results and results[0].get('pmcid'):
            pmcid=results[0]['pmcid']
            try: save(name+'.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmcid+'/fullTextXML')
            except Exception as e: print(name,'full text retrieval',type(e).__name__,str(e))
for name,gid in [('LGSN','51557'),('PLIN4','729359'),('TGM6','343641'),('GRM2','2912'),('Or8b47','258873')]:
    try:
        data=save(name+'_ncbi_gene.json','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&id='+gid)
        print(name,'metadata bytes',len(data))
    except Exception as e: print(name,'metadata error',type(e).__name__,str(e))
    time.sleep(.4)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
