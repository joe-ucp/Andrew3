"""Prespecified descriptive raw-guide robustness; mouse is the biological unit.

This does not reproduce the publication's negative-binomial interaction model.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ['HIVEP2','RNLS','GCN1','MED11','ZBED3','ZNF28']

def cpm(counts):
    totals = counts.sum(axis=0)
    if (totals <= 0).any():
        raise ValueError('Cannot normalize a zero-total sample')
    return counts.divide(totals, axis=1) * 1e6

def contrast(values, treated, controls):
    return values[treated].mean(axis=1) - values[controls].mean(axis=1)

def theta(values, selected):
    n, m = len(values), int(np.sum(selected))
    if not 0 < m < n: raise ValueError('Nonempty candidate and background required')
    r = pd.Series(values).rank(method='average').to_numpy()
    return float((r[selected].sum()-m*(m+1)/2)/(m*(n-m)))

def mouse_eligibility(raw):
    """Current MGI homolog groups plus unique synonym resolution to old library labels."""
    hom = pd.read_csv(raw/'HOM_MouseHumanSequence.rpt',sep='\t',dtype=str).fillna('')
    markers = pd.read_csv(raw/'MRK_List2.rpt',sep='\t',dtype=str).fillna('')
    mouse = hom.loc[hom['NCBI Taxon ID'].eq('10090')]
    human = hom.loc[hom['NCBI Taxon ID'].eq('9606')]
    lib = pd.read_csv(raw/'mouse_geckov2_library_a_2.csv',dtype=str)
    old = set(lib.gene_id)
    valid = set(mouse['Mouse MGI ID'])
    aliases = {}
    for row in markers.loc[markers['MGI Accession ID'].isin(valid)].to_dict('records'):
        for alias in row['Marker Synonyms (pipe-separated)'].split('|'):
            if alias: aliases.setdefault(alias,set()).add(row['Marker Symbol'])
    resolved = set(old)
    alias_used = {}
    current_symbols = set(mouse.Symbol)
    for symbol in old:
        # Exact current symbol takes precedence over historical synonym.
        if symbol not in current_symbols and len(aliases.get(symbol,())) == 1:
            current = next(iter(aliases[symbol])); resolved.add(current)
            alias_used.setdefault(current,[]).append(symbol)
    eligible_rows = mouse[mouse.Symbol.isin(resolved)]
    groups = {}
    for row in eligible_rows.to_dict('records'):
        groups.setdefault(row['DB Class Key'],[]).append(row['Symbol'])
    rows=[]
    for row in human.to_dict('records'):
        m = groups.get(row['DB Class Key'],[])
        rows.append({'gene':row['Symbol'],'human_entrez_id':row['EntrezGene ID'],
                     'mgi_homology_class':row['DB Class Key'],
                     'eligible_mouse_symbols':';'.join(sorted(m)),
                     'library_aliases':';'.join(sorted({a for s in m for a in alias_used.get(s,[])})),
                     'mouse_library_eligible':bool(m), 'eligible_mouse_members':len(m)})
    return pd.DataFrame(rows)

def main():
    raw=ROOT/'data/raw/counts'; processed=ROOT/'data/processed'; results=ROOT/'results'
    processed.mkdir(exist_ok=True); results.mkdir(exist_ok=True)
    pd.read_json(ROOT/'protocol/raw_sensitivity.json',typ='series')  # Require the freeze record.
    counts=pd.read_csv(raw/'GSM6008413_counts-JD_GPP2373_2901507_Sintov_20210224.txt.gz',sep='\t')
    assert len(counts)==77441 and counts['Construct Barcode'].is_unique
    ref=pd.read_csv(raw/'GSM6008413_CP0043_reference_20160120.csv.gz',header=None,names=['barcode','construct','sequence'],dtype=str)
    assert ref.barcode.is_unique
    chk=counts.merge(ref,left_on='Construct Barcode',right_on='barcode',validate='one_to_one')
    assert len(chk)==77441 and (chk['Construct IDs']==chk['construct']).all()
    ann=pd.read_csv(raw/'CP0043_origtarget_20220324.chip',sep='\t',dtype=str)
    ann.columns=['barcode','gene','gene_id']
    strict=pd.read_csv(raw/'CP0043_GRCh38_NCBI_CRISPRko_strict_gene_20240130.chip',sep='\t',dtype=str)
    strict.columns=['barcode','strict_gene','strict_gene_id']
    # Preserve original target identity. Current-score aliases are not guessed.
    mapped=counts.merge(ann,left_on='Construct Barcode',right_on='barcode',validate='one_to_one')
    assert len(mapped)==77441
    strict_audit=mapped[['barcode','gene','gene_id']].merge(strict,on='barcode',how='left')
    strict_audit['same_gene_id']=strict_audit.gene_id.eq(strict_audit.strict_gene_id)
    strict_audit.to_csv(processed/'raw_annotation_audit.csv',index=False)
    metadata=pd.read_csv(raw/'GSM6008413_sample_info.txt.gz',sep='\t')
    control=metadata.loc[metadata['characteristics: sample_type'].eq('Control grafts'),'Sample name'].tolist()
    treated=metadata.loc[metadata['characteristics: sample_type'].eq('hPi grafts'),'Sample name'].tolist()
    assert set(control)=={'Ms #64 Ctrl','Ms #69 Ctrl','Ms #70 Ctrl','Ms #73 Ctrl','Ms #74 Ctrl','Ms #78 Ctrl'}
    assert set(treated)=={'Ms #65','Ms #67','Ms #68','Ms #72','Ms #76','Ms #79'}
    samples=control+treated+['DA446 Pre-Tx','CP0043 Library DNA']
    x=mapped[samples].astype(float)
    assert x.ge(0).all().all() and (x==np.floor(x)).all().all()
    normalized=cpm(x); logs=np.log2(normalized+1)
    sample_qc=pd.DataFrame({'sample':samples,'total_reads':x.sum().to_numpy(),
                            'zero_guides':x.eq(0).sum().to_numpy(),
                            'zero_fraction':x.eq(0).mean().to_numpy()})
    sample_qc['sample_class']=['no_PBMC']*6+['hPi']*6+['pre_transplant','library_DNA']
    sample_qc.to_csv(results/'raw_sample_qc.csv',index=False)
    gs=logs.groupby(mapped.gene).mean()
    cov=normalized.groupby(mapped.gene).mean()
    gene=pd.DataFrame({'guide_count':mapped.groupby('gene').size(),
                        'mean_control_CPM':cov[control].mean(axis=1),
                        'mean_hPi_CPM':cov[treated].mean(axis=1),
                        'mean_pretransplant_CPM':cov['DA446 Pre-Tx'],
                        'mean_library_CPM':cov['CP0043 Library DNA'],
                        'descriptive_log2_CPM_contrast':contrast(gs,treated,control)})
    gene.to_csv(processed/'raw_gene_covariates.csv',index_label='gene')
    gs.to_csv(processed/'raw_gene_sample_log2CPM.csv',index_label='gene')
    selected=mapped.gene.isin(CANDIDATES)
    guides=[]
    for i,row in mapped.loc[selected].iterrows():
        vals={'gene':row.gene,'barcode':row.barcode,'construct_id':row['Construct IDs'],
              'gene_id':row.gene_id,
              'guide_log2_CPM_contrast':float(logs.loc[i,treated].mean()-logs.loc[i,control].mean())}
        for s in samples:
            vals[s+'_count']=int(x.loc[i,s]); vals[s+'_CPM']=float(normalized.loc[i,s])
        guides.append(vals)
    pd.DataFrame(guides).to_csv(results/'raw_candidate_guides.csv',index=False)
    looa=[];loog=[];summary=[]
    scores=pd.read_csv(processed/'allogeneic_scores.csv').set_index('gene')
    for g in CANDIDATES:
        if g not in gs.index:
            summary.append({'gene':g,'status':'unmapped_in_raw_annotation'});continue
        idx=mapped.index[mapped.gene.eq(g)]
        ar=[];gr=[]
        for animal in control+treated:
            t=[s for s in treated if s!=animal]; c=[s for s in control if s!=animal]
            v=float(contrast(gs.loc[[g]],t,c).iloc[0]);ar.append(v)
            looa.append({'gene':g,'omitted_animal':animal,'contrast':v})
        for guide in idx:
            use=idx[idx!=guide]
            v=float(logs.loc[use,treated].mean().mean()-logs.loc[use,control].mean().mean())
            gr.append(v);loog.append({'gene':g,'omitted_barcode':mapped.loc[guide,'barcode'],'contrast':v})
        base=float(gene.loc[g,'descriptive_log2_CPM_contrast'])
        pub=float(scores.loc[g,'published_score'])
        a=strict_audit.loc[strict_audit.gene.eq(g)]
        summary.append({'gene':g,'status':'mapped','guide_count':len(idx),
                        'published_score':pub,'descriptive_log2_CPM_contrast':base,
                        'same_direction_as_published':bool(np.sign(base)==np.sign(pub)),
                        'animal_LOO_min':min(ar),'animal_LOO_max':max(ar),
                        'animal_LOO_all_same_sign':bool(all(np.sign(v)==np.sign(base) for v in ar)),
                        'guide_LOO_min':min(gr),'guide_LOO_max':max(gr),
                        'guide_LOO_all_same_sign':bool(all(np.sign(v)==np.sign(base) for v in gr)),
                        'strict_gene_id_all_concordant':bool(a.same_gene_id.all()),
                        **gene.loc[g].drop(['guide_count','descriptive_log2_CPM_contrast']).to_dict()})
    pd.DataFrame(summary).to_csv(results/'raw_candidate_summary.csv',index=False)
    pd.DataFrame(looa).to_csv(results/'raw_leave_one_animal_out.csv',index=False)
    pd.DataFrame(loog).to_csv(results/'raw_leave_one_guide_out.csv',index=False)
    # Explicitly post-primary exploratory amendment: isolate reference-control choice.
    nt_barcodes=set(strict.loc[strict.strict_gene.str.startswith('NO_SITE_',na=False),'barcode'])
    nt_idx=mapped.index[mapped.barcode.isin(nt_barcodes)]
    assert len(nt_idx)==941
    nt_mean=logs.loc[nt_idx].mean(axis=0)
    nt_delta=float(nt_mean[treated].mean()-nt_mean[control].mean())
    centered=gs.subtract(nt_mean,axis=1)
    nt_rows=[]
    for g in CANDIDATES:
        idx=mapped.index[mapped.gene.eq(g)]
        base=float(contrast(centered.loc[[g]],treated,control).iloc[0])
        ar=[];gr=[]
        for animal in control+treated:
            t=[s for s in treated if s!=animal];c=[s for s in control if s!=animal]
            ar.append(float(contrast(centered.loc[[g]],t,c).iloc[0]))
        for guide in idx:
            use=idx[idx!=guide]
            gr.append(float(logs.loc[use,treated].mean().mean()-logs.loc[use,control].mean().mean()-nt_delta))
        nt_rows.append({'gene':g,'NT_guides':len(nt_idx),'NT_group_difference':nt_delta,
                        'NT_centered_descriptive_contrast':base,'published_score':float(scores.loc[g,'published_score']),
                        'same_direction_as_published':bool(np.sign(base)==np.sign(scores.loc[g,'published_score'])),
                        'animal_LOO_min':min(ar),'animal_LOO_max':max(ar),
                        'guide_LOO_min':min(gr),'guide_LOO_max':max(gr),
                        'analysis_status':'post_primary_exploratory_different_estimand'})
    pd.DataFrame(nt_rows).to_csv(results/'raw_exploratory_NT_centered.csv',index=False)
    pd.DataFrame({'sample':samples,'NT_mean_log2_CPM_plus1':nt_mean.to_numpy()}).to_csv(results/'raw_exploratory_NT_sample_reference.csv',index=False)
    # Prespecified non-causal covariate comparison, against scored genes only.
    universe=scores.join(gene,how='left')
    complete=universe.dropna(subset=['guide_count','mean_control_CPM','mean_pretransplant_CPM']).copy()
    complete['control_quintile']=pd.qcut(complete.mean_control_CPM,5,labels=False,duplicates='drop')
    complete['pretransplant_quintile']=pd.qcut(complete.mean_pretransplant_CPM,5,labels=False,duplicates='drop')
    matched=[]
    for g in CANDIDATES:
        if g not in complete.index:
            matched.append({'gene':g,'status':'unmapped'});continue
        r=complete.loc[g]
        peers=complete.loc[complete.guide_count.eq(r.guide_count)&complete.control_quintile.eq(r.control_quintile)&complete.pretransplant_quintile.eq(r.pretransplant_quintile)&~complete.index.isin(CANDIDATES)]
        v=peers.published_score.to_numpy();s=r.published_score
        percentile=float(np.mean((s>v)+.5*(s==v))) if len(v) else None
        matched.append({'gene':g,'matched_background_n':len(peers),'control_quintile':int(r.control_quintile),
                         'pretransplant_quintile':int(r.pretransplant_quintile),'guide_count':int(r.guide_count),
                         'candidate_vs_matched_probability':percentile,
                         'status':'descriptive_adequate_stratum' if len(peers)>=20 else 'unresolved_small_stratum'})
    pd.DataFrame(matched).to_csv(results/'raw_matched_background.csv',index=False)
    universe.loc[universe.guide_count.isna()].to_csv(results/'raw_scored_genes_without_covariates.csv',index_label='gene')
    # Full-universe mouse library availability, not capitalization as orthology.
    eligibility=mouse_eligibility(raw)
    eligibility.to_csv(processed/'raw_mouse_eligibility.csv',index=False)
    eligible=set(eligibility.loc[eligibility.mouse_library_eligible,'gene'])
    restricted=scores.loc[scores.index.isin(eligible)].copy()
    membership=restricted.index.isin(CANDIDATES)
    coverage=[g for g in CANDIDATES if g in restricted.index]
    mouse_result={'scored_universe_n':len(scores),'eligible_scored_n':len(restricted),
                  'candidate_genes_retained':coverage,'candidate_genes_lost':[g for g in CANDIDATES if g not in coverage],
                  'theta':theta(restricted.published_score.to_numpy(),membership) if coverage else None,
                  'inference':'descriptive restricted-background sensitivity; no replacement of primary test',
                  'boundary':'MGI homology classes can include multiple members; availability is not equivalent to comparable efficacy, expression or discovery survival; no full autoimmune count matrix.'}
    (results/'raw_mouse_library_sensitivity.json').write_text(json.dumps(mouse_result,indent=2))
    complete.to_csv(processed/'raw_scored_universe_covariates.csv',index_label='gene')
    overall={'raw_barcodes':len(counts),'annotation_mapped_barcodes':len(mapped),
             'graft_animals':len(control)+len(treated),'no_PBMC_animals':len(control),'hPi_animals':len(treated),
             'scored_genes_with_raw_covariates':len(complete),'scored_genes_missing_covariates':len(universe)-len(complete),
             'excluded_columns':['Unnamed: 2 (unlabeled zero column)','NTC (negative control; not a graft)'],
             'donor_adjustment':'unavailable: sample table provides no donor assignments',
             'count_model_reproduction':False,
             'warning':'CPM contrast is not the publication KO-by-PBMC interaction relative to non-targeting guides. Differences are diagnostic, not an independent biological replication.'}
    (results/'raw_audit_summary.json').write_text(json.dumps(overall,indent=2))
    print(pd.DataFrame(summary)[['gene','descriptive_log2_CPM_contrast','same_direction_as_published','animal_LOO_min','animal_LOO_max','guide_LOO_min','guide_LOO_max']].to_string(index=False))
    print(json.dumps(mouse_result,indent=2))

if __name__=='__main__': main()
