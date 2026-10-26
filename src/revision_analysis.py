"""Secondary calibration of the primary nominated-set analysis.

Conditional reference distributions concern gene labels and published point
estimates. They do not propagate biological-effect uncertainty.
"""
from pathlib import Path
import gzip
import hashlib
import io
import json
import math
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from src.transfer import rankdata,sample_subsets,superiority
from src.portable_tables import verify_revision_inputs


def theta_from_rank_sums(sums,n,m):
    return (np.asarray(sums)-m*(m+1)/2)/(m*(n-m))


def uniform_theta_draws(scores,m,draws,seed,batch=50_000):
    ranks=rankdata(scores);n=len(ranks)
    rng=np.random.Generator(np.random.PCG64(seed));out=np.empty(draws)
    for start in range(0,draws,batch):
        size=min(batch,draws-start)
        rows=sample_subsets(rng,n,m,size)
        out[start:start+size]=theta_from_rank_sums(ranks[rows].sum(axis=1),n,m)
    return out


def mc_upper_tail(draws,observed):
    exceed=int(np.count_nonzero(draws>=observed));n=len(draws)
    p=(exceed+1)/(n+1)
    return {'draws':n,'exceedances':exceed,'p_upper_tail_plus_one':p,
            'monte_carlo_se':math.sqrt(p*(1-p)/n)}


def no_hit_probability(n,k,m):
    if not(0<=k<=n and 0<=m<=n): raise ValueError('Invalid finite-population counts')
    if n-k<m: return 0.
    return math.prod((n-k-j)/(n-j) for j in range(m))


def stratified_theta(values,strata,selected):
    """Candidate-weighted within-stratum U; no comparison across unmatched strata."""
    total=int(np.sum(selected));result=0.
    for key in sorted(set(strata)):
        group=np.array([x==key for x in strata])
        use=np.flatnonzero(selected[group]);k=len(use)
        if not k: continue
        result+=(k/total)*superiority(values[group],use)
    return float(result)


def matched_randomization(frame,candidates,draws,seed):
    keys=['guide_count','control_quintile','pretransplant_quintile']
    selected=frame.index.isin(candidates);m=int(selected.sum())
    rng=np.random.Generator(np.random.PCG64(seed))
    null=np.zeros(draws);observed=0.;strata=[]
    # Exact tuple keys partition the complete-case universe. Sampling strata
    # independently cannot create cross-stratum collisions.
    for key,group in frame.groupby(keys,sort=True):
        members=np.flatnonzero(group.index.isin(candidates));k=len(members)
        if not k: continue
        values=group.published_score.to_numpy();ranks=rankdata(values);n=len(group)
        if k==n: raise ValueError('A nominated stratum has no background genes')
        t=superiority(values,members);observed+=(k/m)*t
        strata.append({'guide_count':int(key[0]),'control_quintile':int(key[1]),
                       'pretransplant_quintile':int(key[2]),'stratum_size':n,
                       'candidate_count':k,'background_count':n-k,
                       'candidate_genes':';'.join(group.index[members]),
                       'observed_within_stratum_theta':t,'weight':k/m})
        for start in range(0,draws,50_000):
            size=min(50_000,draws-start)
            rows=sample_subsets(rng,n,k,size)
            null[start:start+size]+=(k/m)*theta_from_rank_sums(ranks[rows].sum(axis=1),n,k)
    return float(observed),null,pd.DataFrame(strata)


def weighted_rank_tuples(rng,n,alphas,size):
    """Independent Beta(a_i,1) rank quantiles, conditional on distinct genes.

    Position i selects zero-based finite rank j with probability
    ((j+1)/N)**a_i - (j/N)**a_i before distinctness conditioning.
    Resampling an entire colliding row preserves that conditional joint law.
    """
    alphas=np.asarray(alphas,dtype=float)
    if n<len(alphas) or (alphas<=0).any(): raise ValueError('Invalid rank-alternative parameters')
    def draw(k):
        u=rng.power(alphas,size=(k,len(alphas)))
        return np.minimum((u*n).astype(np.int64),n-1)
    rows=draw(size);rejected=0
    while True:
        bad=np.any(np.diff(np.sort(rows,axis=1),axis=1)==0,axis=1)
        count=int(bad.sum())
        if not count: return rows,rejected
        rejected+=count;rows[bad]=draw(count)


def save_draws(path,draws):
    data=io.BytesIO();np.save(data,draws,allow_pickle=False)
    path.write_bytes(gzip.compress(data.getvalue(),mtime=0))


def save_histogram(path,draws):
    counts,edges=np.histogram(draws,bins=np.linspace(0,1,201))
    pd.DataFrame({'bin_left':edges[:-1],'bin_right':edges[1:],'count':counts,
                  'density':counts/(len(draws)*np.diff(edges))}).to_csv(path,sep='\t',index=False,float_format='%.17g')


def main():
    plan=json.loads((ROOT/'protocol/revision_analysis_plan.json').read_text(encoding='utf-8'))
    verify_revision_inputs(ROOT,plan['input_sha256'])
    out=ROOT/'results';out.mkdir(exist_ok=True)
    d=pd.read_csv(ROOT/'data/processed/allogeneic_scores.csv').set_index('gene')
    candidates=plan['candidates'];scores=d.published_score.to_numpy()
    idx=np.flatnonzero(d.index.isin(candidates));n=len(d);m=len(idx)
    assert n==plan['primary_null']['n'] and m==plan['primary_null']['m']
    null=uniform_theta_draws(scores,m,plan['primary_null']['draws'],plan['primary_null']['seed'])
    observed=superiority(scores,idx)
    quantiles=np.quantile(null,plan['primary_null']['quantiles'],method='higher')
    quantile_map=dict(zip(plan['primary_null']['quantiles'],quantiles))
    threshold=float(quantile_map[.95])
    qtable=pd.DataFrame({'quantile':plan['primary_null']['quantiles'],'theta':quantiles,
                        'meaning':'empirical competitive-null reference, not biological confidence limit'})
    qtable.to_csv(out/'revision_primary_null_quantiles.tsv',sep='\t',index=False,float_format='%.17g')
    save_histogram(out/'revision_primary_null_histogram.tsv',null)
    save_draws(out/'revision_primary_null.npy.gz',null)
    primary={'universe_n':n,'candidates_m':m,'observed_theta':observed,
             **mc_upper_tail(null,observed),'seed':plan['primary_null']['seed'],
             'null_mean':float(null.mean()),'null_sd':float(null.std(ddof=1)),
             'central_95_lower':float(quantile_map[.025]),'central_95_upper':float(quantile_map[.975]),
             'upper_5pct_reference_boundary':threshold,'rejection_rule':'theta > empirical .95 quantile',
             'reference_null_rejection_fraction':float(np.mean(null>threshold)),
             'threshold_tie_count':int(np.sum(null==threshold)),
             'interpretation':'Conditional published-point-estimate/gene-label reference; not a biological-effect interval'}
    k=int(((d.published_score>0)&(d.reported_q<=.05)).sum())
    nohit={'universe_n':n,'positive_genomewide_supported_genes':k,'candidate_set_size':m,
           'observed_hits':int(((d.loc[candidates].published_score>0)&(d.loc[candidates].reported_q<=.05)).sum()),
           'uniform_random_set_probability_zero_hits':no_hit_probability(n,k,m),
           'expected_hit_count':m*k/n,'formula':'C(N-K,m)/C(N,m)',
           'interpretation':'Calibration of discovery-threshold count; not an independent negative biological test'}
    common=pd.read_csv(ROOT/'data/processed/raw_scored_universe_covariates.csv').set_index('gene')
    # Original strata are reused verbatim, and checked against their frozen rule.
    assert common.index.is_unique and set(candidates)<=set(common.index)
    for cov,column in [('mean_control_CPM','control_quintile'),('mean_pretransplant_CPM','pretransplant_quintile')]:
        recomputed=pd.qcut(common[cov],5,labels=False,duplicates='drop')
        assert np.array_equal(recomputed.to_numpy(),common[column].to_numpy())
    matched_observed,matched_null,strata=matched_randomization(common,candidates,
        plan['joint_matching']['draws'],plan['joint_matching']['seed'])
    strata.to_csv(out/'revision_matched_strata.tsv',sep='\t',index=False,float_format='%.17g')
    save_histogram(out/'revision_matched_null_histogram.tsv',matched_null)
    save_draws(out/'revision_matched_null.npy.gz',matched_null)
    matched={'common_universe_n':len(common),'universe_missing_covariates':len(d)-len(common),
             'candidate_count':m,'n_occupied_strata':len(strata),'observed_T':matched_observed,
             **mc_upper_tail(matched_null,matched_observed),'seed':plan['joint_matching']['seed'],
             'null_mean':float(matched_null.mean()),'null_sd':float(matched_null.std(ddof=1)),
             'central_95_lower':float(np.quantile(matched_null,.025,method='higher')),
             'central_95_upper':float(np.quantile(matched_null,.975,method='higher')),
             'reference':'preserved nominee count per disjoint stratum, full stratum includes nominees',
             'scope':'Secondary conditional label reference; biological gene correlations unmodeled'}
    samples=pd.read_csv(ROOT/'results/raw_exploratory_NT_sample_reference.csv').set_index('sample')
    control=[s for s in samples.index if s.startswith('Ms #') and s.endswith('Ctrl')]
    treated=[s for s in samples.index if s.startswith('Ms #') and not s.endswith('Ctrl')]
    assert len(control)==len(treated)==6
    nt=float(samples.loc[treated,'NT_mean_log2_CPM_plus1'].mean()-samples.loc[control,'NT_mean_log2_CPM_plus1'].mean())
    cidx=np.flatnonzero(common.index.isin(candidates))
    published=d.loc[common.index,'published_score'].to_numpy()
    raw=common.descriptive_log2_CPM_contrast.to_numpy();centered=raw-nt
    common_rows=[]
    for name,values in [('published_score',published),('raw_log2_CPM_contrast',raw),('raw_NT_centered',centered)]:
        common_rows.append({'score_definition':name,'universe_n':len(common),'candidate_count':m,
                            'theta':superiority(values,cidx),'positive_candidate_scores':int((values[cidx]>0).sum()),
                            'scope':'exploratory point-estimate comparison; not exact source-model reconstruction'})
    common_table=pd.DataFrame(common_rows)
    common_table.to_csv(out/'revision_common_universe.tsv',sep='\t',index=False,float_format='%.17g')
    shift={'shared_NT_constant':nt,'rank_vectors_identical':bool(np.array_equal(rankdata(raw),rankdata(centered))),
           'theta_identical':bool(common_rows[1]['theta']==common_rows[2]['theta']),
           'identity':'(d_g-c)>(d_h-c) iff d_g>d_h; equalities are also preserved',
           'consequence':'Common centering changes the reference zero/sign, not the ordering or rank enrichment in a fixed universe'}
    assert shift['rank_vectors_identical'] and shift['theta_identical']
    # Frozen stylized label alternatives; no observed-effect-derived parameter.
    ordered_ranks=np.sort(rankdata(scores));alt=[]
    count=plan['stylized_alternatives']['draws_per_scenario']
    for scenario in plan['stylized_alternatives']['scenarios']:
        alphas=np.array(scenario['alphas'],dtype=float)
        rng=np.random.Generator(np.random.PCG64(scenario['seed']))
        vals=np.empty(count);rejected=0
        for start in range(0,count,50_000):
            size=min(50_000,count-start)
            rows,r=weighted_rank_tuples(rng,n,alphas,size);rejected+=r
            vals[start:start+size]=theta_from_rank_sums(ordered_ranks[rows].sum(axis=1),n,m)
        detections=int((vals>threshold).sum());rate=detections/count
        alt.append({'scenario':scenario['name'],'alpha_vector':';'.join(map(str,alphas)),
                    'nominal_mean_rank_quantile_before_conditioning':float(np.mean(alphas/(alphas+1))),
                    'simulated_mean_theta':float(vals.mean()),'mean_theta_monte_carlo_se':float(vals.std(ddof=1)/math.sqrt(count)),
                    'draws':count,'seed':scenario['seed'],'rejected_duplicate_tuples':rejected,
                    'reference_rejection_boundary':threshold,'rejections':detections,
                    'rejection_fraction':rate,'rejection_fraction_monte_carlo_se':math.sqrt(rate*(1-rate)/count),
                    'interpretation':'stylized finite-rank alternative; not observed power or experimental-effect uncertainty'})
    pd.DataFrame(alt).to_csv(out/'revision_alternative_sensitivity.tsv',sep='\t',index=False,float_format='%.17g')
    summary={'revision_plan_sha256':hashlib.sha256((ROOT/'protocol/revision_analysis_plan.json').read_bytes()).hexdigest(),
             'primary_null_calibration':primary,'no_hit_calibration':nohit,
             'secondary_joint_matched_test':matched,'common_universe_constant_shift':shift,
             'experimental_uncertainty':'Not estimated: published vector supplies point effects and adjusted q, not nominal p-values, covariance or complete model uncertainty. No inverse-q reconstruction attempted.',
             'calculation_scope':'All numerical results computed from the retained scientific inputs.'}
    (out/'revision_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    print(common_table.to_string(index=False))
    print(pd.DataFrame(alt)[['scenario','simulated_mean_theta','rejection_fraction','rejection_fraction_monte_carlo_se']].to_string(index=False))


if __name__=='__main__': main()
