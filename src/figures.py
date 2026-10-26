"""Publication figures from computed results and explicit source curation."""
from pathlib import Path
import json
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT=Path(__file__).resolve().parents[1]
NAVY='#193C59'; TEAL='#247C88'; GOLD='#B36C31'; GRAY='#6E7E89'

def save(fig,name):
    out=ROOT/'manuscript/figures';out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/(name+'.png'),dpi=220,facecolor='white')
    fig.savefig(out/(name+'.svg'),metadata={'Date':None});plt.close(fig)

def box(ax,xy,w,h,title,body,color='#EEF3F6',fontsize=12):
    x,y=xy
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.016',facecolor=color,edgecolor='#AFBEC8',linewidth=1))
    ax.text(x+.018,y+h-.032,title,weight='bold',fontsize=13,va='top',color=NAVY)
    ax.text(x+.018,y+h-.095,body,fontsize=fontsize,va='top',linespacing=1.5)

def eligibility():
    mapping=pd.read_csv(ROOT/'candidate_mapping.tsv',sep='\t')
    assert len(mapping)==11 and int(mapping.included.sum())==6
    fig,ax=plt.subplots(figsize=(10,6));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    box(ax,(.02,.69),.45,.28,'A  Autoimmune discovery | Cai 2020','Mouse NIT-1 cells • subcutaneous grafts\nNOD.scid hosts + diabetic NOD splenocytes\nPooled CRISPR screen; retrieval at 8 weeks')
    box(ax,(.53,.69),.45,.28,'B  Allogeneic evaluation | Sintov 2022','Human SC-islets • kidney-capsule grafts\nMHC-null NSG hosts + unmatched human PBMCs\nPublished gene scores; retrieval at 10 weeks')
    box(ax,(.02,.405),.96,.20,'11 original discovery labels','Hivep2 • Rnls • Gcn1l1 • Med11 • Zbed3 • Gm3604\nLgsn • Plin4 • Tgm6 • Grm2 • Olfr909',fontsize=12.5)
    ax.annotate('',xy=(.24,.612),xytext=(.24,.68),arrowprops=dict(arrowstyle='-|>',color=GRAY))
    box(ax,(.02,.09),.43,.235,'6 scored human counterparts','HIVEP2 • RNLS • GCN1\nMED11 • ZBED3 • ZNF28*',color='#E6F1F1',fontsize=13)
    box(ax,(.53,.09),.45,.235,'5 outside quantitative comparison','4 mapped but unscored: LGSN, PLIN4,\nTGM6, GRM2\n1 uncommitted mapping: Olfr909',fontsize=11.5)
    for x in [.24,.75]: ax.annotate('',xy=(x,.333),xytext=(x,.398),arrowprops=dict(arrowstyle='-|>',color=GRAY))
    ax.text(.02,.025,'All 11 retained in the functional map.  *Gm3604/ZNF28 mapping omitted in a separate sensitivity.',fontsize=11)
    fig.subplots_adjust(left=.02,right=.98,top=.99,bottom=.01);save(fig,'eligibility')

def quantitative():
    c=pd.read_csv(ROOT/'results/candidate_scores.tsv',sep='\t');d=pd.read_csv(ROOT/'data/processed/allogeneic_scores.csv')
    t=pd.read_csv(ROOT/'results/transfer_summary.tsv',sep='\t').set_index('analysis')
    summary=json.loads((ROOT/'results/revision_summary.json').read_text())
    common=pd.read_csv(ROOT/'results/revision_common_universe.tsv',sep='\t').set_index('score_definition')
    null=summary['primary_null_calibration'];hist=pd.read_csv(ROOT/'results/revision_primary_null_histogram.tsv',sep='\t')
    fig=plt.figure(figsize=(10,7.1));outer=fig.add_gridspec(2,1,height_ratios=[1.4,1],hspace=.95)
    top=outer[0].subgridspec(1,2,wspace=.40)
    ax=fig.add_subplot(top[0]);scores=np.sort(d.published_score.to_numpy())
    for i,(_,row) in enumerate(c.iterrows()):
        pct=100*((scores<row.published_score).sum()+.5*(scores==row.published_score).sum())/len(scores)
        ax.plot([0,100],[i,i],color='#E5EBEF',lw=6,zorder=1)
        ax.scatter(pct,i,s=48,color=TEAL if row.published_score>0 else GOLD,zorder=4)
        ax.text(pct+3,i,f'{pct:.1f}',va='center',fontsize=10)
    ax.axvline(50,color='#BBB',lw=.8,ls=':');ax.set_xlim(0,104);ax.set_xticks([0,25,50,75,100])
    ax.set_yticks(range(len(c)),c.gene);ax.invert_yaxis();ax.set_xlabel('Percentile in all 11,689 scores')
    ax.set_title('A  Positions in the reported ordering',loc='left',weight='bold',fontsize=12,pad=14)
    ax=fig.add_subplot(top[1]);mid=(hist.bin_left+hist.bin_right)/2
    ax.bar(mid,hist.density,width=hist.bin_right-hist.bin_left,color='#C5D2DA',edgecolor='none')
    crit=null['upper_5pct_reference_boundary'];obs=null['observed_theta']
    ax.axvspan(crit,1,color='#EAD6C4',alpha=.7)
    ax.axvline(obs,color=NAVY,lw=2,label=f'Observed {obs:.3f}')
    ax.axvline(.5,color=GRAY,ls=':',lw=1.6,label='Neutral reference 0.5')
    ax.axvline(crit,color=GOLD,ls='--',lw=1.6,label=f'Upper 5% boundary {crit:.3f}')
    ax.set(xlim=(0,1),xlabel='Six-gene statistic, θ',ylabel='Null density')
    ax.set_title('B  Competitive null calibration',loc='left',weight='bold',fontsize=12,pad=14)
    ax.legend(fontsize=9.5,frameon=False,loc='upper left',bbox_to_anchor=(-.02,-.27))
    bottom=outer[1].subgridspec(1,3,wspace=.75)
    groups=[
        ('C1  Candidate composition',['Six genes','No ZNF28'],[obs,t.loc['exclude_ZNF28','theta']],[NAVY,TEAL],'11,689-gene universe','Superiority, θ'),
        ('C2  Score definition',['Published','Raw / NT-centered'],[common.loc['published_score','theta'],common.loc['raw_log2_CPM_contrast','theta']],[NAVY,GOLD],'Shared 11,474-gene universe','Superiority, θ'),
        ('C3  Matched reference',['Matched set'],[summary['secondary_joint_matched_test']['observed_T']],[GRAY],'Coverage/abundance strata','Matched statistic, T')]
    for j,(title,labels,values,colors,note,xlab) in enumerate(groups):
        ax=fig.add_subplot(bottom[j]);ys=np.arange(len(values))
        ax.scatter(values,ys,s=55,c=colors,zorder=3)
        for i,value in enumerate(values):ax.text(value+.008,i,f'{value:.3f}',va='center',fontsize=10)
        ax.axvline(.5,color='#AAA',ls=':',lw=1);ax.set(xlim=(.32,.57),ylim=(len(values)-.5,-.8),xlabel=xlab)
        ax.set_yticks(ys,labels,fontsize=9.5);ax.set_xticks([.35,.45,.55]);ax.set_title(title,loc='left',fontsize=11,weight='bold',pad=19)
        ax.text(0,1.04,note,transform=ax.transAxes,fontsize=8.5)
        ax.set_facecolor(['#F5F8FA','#F7F5F2','#F4F5F5'][j])
    fig.text(.09,.025,'Panels C1–C3 address different questions; they are not interchangeable estimates of one biological effect.',fontsize=10)
    fig.subplots_adjust(left=.10,right=.97,top=.94,bottom=.13);save(fig,'transfer')

def functional():
    f=pd.read_csv(ROOT/'evidence/functional_matrix.csv',keep_default_na=False);assert len(f)==11
    measured=set(pd.read_csv(ROOT/'evidence/baseline_experiments.tsv',sep='\t').gene)
    fig,ax=plt.subplots(figsize=(10,6.2));ax.axis('off');rows=[]
    for _,r in f.iterrows():
        label=r.gene
        if r.original_mouse_label in ['Gcn1l1','Gm3604','Olfr909']: label+='\n('+r.original_mouse_label+')'
        baseline='Measured\nSee Figure 4' if r.gene in measured else 'NL | search log'
        eligible=str(r.quantitative_eligible).lower()=='true'
        allo='NM | S Fig1' if eligible else 'NL | search log'
        protection=r.get('figure_protection_text','') or {'HIVEP2':'No durable mouse allo\nK Fig8','RNLS':'Auto validated: C Fig2\nMouse allo: C ED3; K Fig8'}.get(r.gene,'Pooled nomination\nC Fig1')
        rows.append([label,baseline,'NM | C Fig1'+ ('/2' if r.gene=='RNLS' else ''),allo,protection])
    headers=['Target','Baseline\nsecretion','Post-autoimmune\nsecretion','Post-allogeneic\nsecretion','Separate protection\nevidence']
    tab=ax.table(cellText=rows,colLabels=headers,colWidths=[.13,.22,.18,.18,.29],cellLoc='center',loc='center')
    tab.auto_set_font_size(False);tab.set_fontsize(10.5);tab.scale(1,2.05)
    for j in range(5):tab[2,j].set_height(tab[2,j].get_height()*1.5)
    for (i,j),cell in tab.get_celld().items():
        cell.set_edgecolor('#CFD7DC');cell.set_linewidth(.65)
        if i==0:cell.set_facecolor('#E4ECF1');cell.set_text_props(weight='bold',fontsize=11)
        else:
            value=rows[i-1][j]
            cell.set_facecolor('#E3F0F1' if value.startswith('Measured') else '#F4F5F6' if value.startswith('NL') else '#E7EDF4' if value.startswith('NM') else 'white')
        if i==7:cell.set_linewidth(1.4)
    ax.set_title('Measurement coverage within the specified eleven-target search',loc='left',weight='bold',fontsize=13,pad=20)
    fig.text(.025,.07,'NM: not measured in an extracted comparison  •  NL: no qualifying comparison located\nBaseline observations: Figure 4. Missingness is not failure. Bode 2023 full text and supplements assessed.',fontsize=10.5,linespacing=1.5)
    fig.text(.025,.015,'C: Cai 2020  •  S: Sintov 2022  •  K: Karaoglu 2026  •  B: Bode 2024  •  ED: Extended Data\nFirst six rows: quantitative cohort. Remaining five: literature coverage. Xenografts kept separate.',fontsize=10,linespacing=1.4)
    fig.subplots_adjust(left=.025,right=.975,top=.88,bottom=.21);save(fig,'evidence_matrix')

def baseline():
    # Figure labels are curated with exact comparison IDs, not inferred from a gene-level grade.
    f=pd.read_csv(ROOT/'evidence/baseline_experiments.tsv',sep='\t',keep_default_na=False)
    columns=['figure_identity','figure_response','figure_output','figure_comparison']
    if not set(columns).issubset(f.columns):
        raise ValueError('Baseline curation must supply the four source-resolved figure labels')
    fig,ax=plt.subplots(figsize=(10,6.9));ax.axis('off')
    wrapped=[['\n'.join(textwrap.fill(line,width=width,break_long_words=False,break_on_hyphens=False) for line in str(value).splitlines()) for value,width in zip(row,[30,31,19,31])] for row in f[columns].values.tolist()]
    table=ax.table(cellText=wrapped,colLabels=['Gene / model / study','Observed glucose response\nand assay control','Output reporting','Comparison with reference'],
                   colWidths=[.26,.28,.18,.28],cellLoc='center',loc='center')
    table.auto_set_font_size(False);table.set_fontsize(10.5);table.scale(1,3.75)
    for (i,j),cell in table.get_celld().items():
        cell.set_edgecolor('#CBD5DD');cell.set_linewidth(.7)
        if i==0:cell.set_facecolor('#E4ECF1');cell.set_text_props(weight='bold')
        else:cell.set_facecolor('white' if i%2 else '#F5F8F9')
    ax.set_title('Baseline secretion findings by gene, model and study',loc='left',fontsize=13,weight='bold',pad=20)
    fig.text(.025,.065,'All rows were measured before or without immune challenge. Observed baseline response is distinct\nfrom comparative equivalence and from post-immune retention. No bar-height values were digitized.',fontsize=10.5,linespacing=1.5)
    fig.text(.025,.015,'Exact record identifiers, assay definitions and source locators: evidence/baseline_experiments.tsv.',fontsize=10)
    fig.subplots_adjust(left=.025,right=.975,top=.87,bottom=.17);save(fig,'baseline_experiments')

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'andrew3-scientific-figures'})
    eligibility();quantitative();functional();baseline()
    print('Wrote eligibility, percentile/calibration and source-specific function figures')

if __name__=='__main__': main()
