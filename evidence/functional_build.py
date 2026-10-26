"""Write the manually adjudicated comparison-level evidence table.

This script serializes source-linked curation; it does not infer unreported assays.
Run from the project root. Numeric values are only transcribed from source text.
"""
from pathlib import Path
import csv
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
N = 'not measured/reported for this comparison'
fields = ['record_id','gene','primary_candidate','study','modification','editing_verification',
          'species','cell_model','challenge_class','effector_source','matched_control',
          'time_relative_to_challenge','biological_n','technical_replication',
          'basal_glucose','basal_insulin','stimulated_glucose','stimulated_insulin',
          'absolute_stimulated_output','stimulation_index','normalization_denominator',
          'uncertainty','control_responsiveness','function_state','adverse_function_flag',
          'immune_validation','protection_result','retention_margin','retention_estimate_interval',
          'source_url','source_locator','source_file','limitations','joint_gate']
rows=[]
urls={'Cai2020':'https://doi.org/10.1038/s42255-020-0254-1',
      'Sintov2022':'https://doi.org/10.1016/j.stemcr.2022.08.002',
      'Karaoglu2026':'https://doi.org/10.3389/fimmu.2026.1759835',
      'MacDonald2025':'https://doi.org/10.1016/j.molmet.2025.102115',
      'Bode2023':'https://doi.org/10.2337/db23-0030',
      'Bode2024':'https://doi.org/10.3389/fimmu.2024.1403752'}
files={'Cai2020':'cai2020','Sintov2022':'sintov2022','Karaoglu2026':'karaoglu2026','MacDonald2025':'rnls2025',
       'Bode2023':'revision3_citations/Bode2023_efetch','Bode2024':'revision3_citations/PMC11225407'}
def add(gene,study,locator,**kw):
    row={k:N for k in fields}
    row.update(record_id=f'FE{len(rows)+1:03d}',gene=gene,primary_candidate=str(gene in candidates).lower(),
               study=study,source_url=urls[study],source_locator=locator,
               source_file=f'data/raw/literature/{files[study]}.xml',
               function_state='F0',adverse_function_flag='not assessed',
               retention_margin='none prespecified in retrieved report',
               retention_estimate_interval='not estimable from extracted data',
               joint_gate='not established')
    row.update(kw);rows.append(row)

candidates=['HIVEP2','RNLS','GCN1','MED11','ZBED3','ZNF28']
mouse={'HIVEP2':'Hivep2','RNLS':'Rnls','GCN1':'Gcn1l1','MED11':'Med11','ZBED3':'Zbed3','ZNF28':'Gm3604'}
for gene in candidates:
    add(gene,'Cai2020','Figure 1; RESULTS: CRISPR screen for protective mutations',
        modification=f'pooled GeCKO-A loss-of-function guide targeting {mouse[gene]}',
        editing_verification='library perturbation; no candidate-specific editing verification in Figure 1',
        species='mouse',cell_model='NIT-1 insulinoma graft in NOD.scid',challenge_class='autoimmune',
        effector_source='splenocytes from diabetic NOD mice',matched_control='graft without diabetogenic splenocytes',
        time_relative_to_challenge='8 weeks after transplantation/co-injection',
        biological_n='not specified in Figure 1 caption; no invented count',
        technical_replication='one detected guide per reported discovery target; not biological replication',
        immune_validation='discovery screen only in this row',protection_result='surviving-guide selection',
        limitations='Guide persistence is not secretion or independent edit validation. Gm3604-to-ZNF28 identity is mapping-sensitive.' if gene=='ZNF28' else 'Guide persistence is not secretion or independent edit validation.')
    add(gene,'Sintov2022','Figure 1G-H; Data S1 sheet Fig3 - in vivo CRISPR screen',
        modification='pooled Brunello CRISPR loss of function',editing_verification='pooled guide counts, not separately verified candidate clone',
        species='human',cell_model='HUES8-derived SC-islet grafts in MHC-null NSG mice',challenge_class='allogeneic',
        effector_source='unmatched healthy human PBMCs (HLA-A2 negative)',matched_control='no-PBMC grafts',
        time_relative_to_challenge='10 weeks after PBMC injection',biological_n='6 graft mice per condition',
        technical_replication='up to four targeting guides per gene; guides are not mice',
        immune_validation='screen only',protection_result='see computed candidate score table; no functional inference',
        limitations='Mixed endocrine population; donor assignments to graft IDs not recovered. Bulk graft insulin cannot be assigned to one pooled edit.')

add('RNLS','Cai2020','Figure 2a-d; Extended Data Figure 2a-b',
    modification='Rnls-mutant NIT-1 line (partial loss of function)',editing_verification='T7E1 and sequencing; ~75% cells predicted to carry at least one deleterious allele',
    species='mouse',cell_model='NIT-1 paired flank grafts in NOD.scid',challenge_class='autoimmune',
    effector_source='diabetic NOD splenocytes',matched_control='non-targeting guide on opposite flank; no-splenocyte recipients',
    time_relative_to_challenge='days 0-18; extended follow-up to day 57',biological_n='5 with splenocytes; 3 without',
    technical_replication='paired grafts within mouse',immune_validation='separately engineered cells',
    protection_result='improved survival in autoimmune setting',uncertainty='mean +/- SEM; source p=0.0357 at days 14/18',
    limitations='Survival via luminescence; no paired basal/stimulated secretion after autoimmune exposure.')
add('RNLS','Cai2020','Figure 2e-h',
    modification='Rnls-mutant NIT-1 line',editing_verification='as FE013',species='mouse',
    cell_model='NIT-1 grafts in overtly diabetic NOD mice',challenge_class='autoimmune',
    effector_source='ongoing endogenous NOD autoimmunity',matched_control='paired non-targeting graft',
    time_relative_to_challenge='days 0-18',biological_n='5 mice',technical_replication='paired grafts',
    immune_validation='separately engineered cells',protection_result='improved survival',
    uncertainty='mean +/- SEM; source p=0.012 for remaining luminescence',
    limitations='No post-challenge basal/stimulated secretion.')
add('RNLS','Cai2020','Figure 2i-l',modification='RIP-Cas9 and Rnls guide in dispersed primary islets',
    editing_verification='targeted lentiviral perturbation; no per-graft edited fraction recovered',species='mouse',
    cell_model='primary NOD.scid islets under paired kidney capsules of NOD.scid recipients',challenge_class='autoimmune',
    effector_source='diabetogenic NOD splenocytes injected two weeks after grafting',matched_control='paired non-targeting islets; no-splenocyte cohort',
    time_relative_to_challenge='grafts assessed day 39; splenocytes at week 2',biological_n='6 autoimmune; 5 non-autoimmune mice',
    technical_replication='paired kidney grafts',immune_validation='primary-islet edit validation',
    protection_result='retained insulin-expression/graft-survival evidence',normalization_denominator='insulin mRNA relative to glucagon mRNA',
    uncertainty='mean +/- SEM; p=0.0312 and p=0.0173 for specified source comparisons',
    limitations='Insulin mRNA is a proxy, not secreted insulin; do not join this row to baseline GSIS as post-immune function.')
add('RNLS','Cai2020','Extended Data Figure 3b-c',modification='Rnls-mutant NIT-1 line',editing_verification='as FE013',
    species='mouse',cell_model='NIT-1 paired flank grafts in C57BL/6',challenge_class='allogeneic',
    effector_source='MHC-mismatched immunocompetent C57BL/6 host',matched_control='paired non-targeting graft',
    time_relative_to_challenge='days 0, 4, 7',biological_n='3 mice',technical_replication='paired grafts',
    immune_validation='separately engineered cells',protection_result='both control and mutant destroyed within a week',
    uncertainty='mean +/- SEM; no equivalence/noninferiority interval',limitations='Known negative context; not a new discovery in this analysis.')
add('RNLS','Cai2020','Extended Data Figure 2c; Methods: GSIS and insulin ELISA',
    modification='Rnls-targeting lentivirus plus RIP-Cas9 in primary islets',editing_verification='delivery documented; assay-specific editing proportion not reported',
    species='mouse',cell_model='dispersed primary CD1 islets',challenge_class='none',effector_source='none',matched_control='non-targeting guide',
    time_relative_to_challenge='72 hours after transduction; no immune challenge',biological_n='donor count not reported; approximately 1700 islets purified',
    technical_replication='5 technical replicates per condition/genotype',basal_glucose='2.8 mM for 1 hour',stimulated_glucose='16.8 mM for 1 hour; then 30 mM KCl',
    basal_insulin='shown graphically; exact values not extracted',stimulated_insulin='shown graphically; exact values not extracted',
    absolute_stimulated_output='not numerically extracted; DNA-normalized ELISA plot',stimulation_index='not reported numerically',
    normalization_denominator='islet genomic DNA content',uncertainty='mean +/- SEM; no significant edited/control difference; no retention interval',
    control_responsiveness='source states dispersion reduced overall responsiveness',function_state='F1',adverse_function_flag='no significant edited/control difference; not noninferiority',
    immune_validation='not an immune assay',protection_result='not applicable',limitations='Technical n and diminished control response preclude broad retention inference.')
add('RNLS','Cai2020','Extended Data Figure 10a-c',modification='clonal biallelic 112 bp RNLS deletion, clone 1',
    editing_verification='PCR and sequencing',species='human',cell_model='iPSC-derived SC-beta cells',challenge_class='none',effector_source='none',
    matched_control='isogenic WT clone',time_relative_to_challenge='baseline differentiated cells; no immune challenge',
    biological_n='caption: four independent SC-beta batches',technical_replication='each batch in triplicate; caption also calls n=4 technical replicates (ambiguity retained)',
    basal_glucose='2.8 mM',stimulated_glucose='20 mM; 30 mM KCl positive stimulus',basal_insulin='shown graphically; not numerically extracted',
    stimulated_insulin='shown graphically; not numerically extracted',absolute_stimulated_output='(microIU/mL)/1000 beta cells shown graphically; exact values not digitized',
    stimulation_index='not numerically extracted',normalization_denominator='1000 beta cells (Extended Data Figure 10c ordinate)',
    uncertainty='mean +/- SEM; no prespecified retention margin/interval',control_responsiveness='WT and RNLS KO show increased insulin at 20 versus 2.8 mM glucose (ED10c visual observation)',
    function_state='F1',adverse_function_flag='no impairment reported in source baseline assay',immune_validation='none in this assay',
    protection_result='not applicable',limitations='ER-stress apoptosis experiment (Figure 6) is separate and is not an allogeneic/autoimmune challenge.')
add('RNLS','MacDonald2025','Extended Data Figure 7D; Methods: Static GSIS; Results: human SC-beta metabolism',
    modification='RNLS clonal knockout; pharmacological arms kept distinct in source',editing_verification='CRISPR clonal line generated as Cai 2020',
    species='human',cell_model='planar iPSC-derived SC-beta cells, stage 6 day 14',challenge_class='none',effector_source='none',
    matched_control='WT SC-beta cells',time_relative_to_challenge='baseline GSIS; no immune exposure',
    biological_n='not clearly declared for panel 7D',technical_replication='3 plotted observations per group; biological independence not assumed',
    basal_glucose='2 mM for 1 hour',stimulated_glucose='20 mM for 1 hour',basal_insulin='plot only; not numerically extracted',
    stimulated_insulin='plot only; not numerically extracted',absolute_stimulated_output='microIU per 100000 cells; values not digitized',
    stimulation_index='not reported numerically; glucose response absent in both WT and KO',normalization_denominator='end-of-assay cell count',
    uncertainty='mean +/- SEM; no retention interval',control_responsiveness='absent GSIS in WT and KO; KCl response present',function_state='F1',
    adverse_function_flag='assay unresponsive; not evidence of retained responsiveness',immune_validation='none in this assay',protection_result='not applicable',
    limitations='Stress-survival and metabolic findings cannot replace glucose response. Additional 2025 follow-up discovered during bounded search.')

for gene in ['HIVEP2','RNLS']:
    for species,model,n,panel in [('mouse','beta-TC-6 spheroids',5,'6A'),('human','EndoC-betaH1 spheroids',3,'6B')]:
        add(gene,'Karaoglu2026','Figure '+panel+'; Methods: GSIS; Figure 3 editing',
            modification=f'CRISPR-Cas9 {gene} guide; antibiotic-selected population',editing_verification='T7E1 and TIDE; protein-level knockout not demonstrated',
            species=species,cell_model=model,challenge_class='none',effector_source='none',matched_control='non-targeting guide spheroids',
            time_relative_to_challenge='pretransplant baseline GSIS',biological_n=f'{n} independent experiments',
            technical_replication='30 spheroids per experiment; not 30 biological replicates',
            basal_glucose='mouse-specific method 1.67 mM; general method 2.8 mM (inconsistent)' if species=='mouse' else '2.8 mM',
            stimulated_glucose='16.7 mM' if species=='mouse' else 'human-specific method 28 mM; general method 16.7 mM (inconsistent)',
            basal_insulin='absolute values not reported in Figure 6',stimulated_insulin='absolute values not reported in Figure 6',
            absolute_stimulated_output='unavailable; only stimulation index shown',stimulation_index='graph only; no raw values extracted',
            normalization_denominator='high-glucose insulin / low-glucose insulin within assay',uncertainty='mean +/- SD; RNLS mouse p<0.05; other candidate contrasts nonsignificant',
            control_responsiveness='NT SI near 1; control unresponsive' if species=='mouse' else 'NT shows glucose stimulation in Figure 6B',
            function_state='F1',adverse_function_flag='reduced SI in mouse RNLS; interpret against unresponsive control' if gene=='RNLS' and species=='mouse' else 'no significant difference; does not establish noninferiority',
            immune_validation='not an immune assay',protection_result='not applicable',limitations='Baseline assay; absolute output/retention interval missing; glucose dose descriptions conflict.')
    for species,model,challenge,locator in [('mouse','beta-TC-6 spheroids in CD1','allogeneic','Figure 8'),('human','EndoC-betaH1 spheroids in CD1','xenogeneic','Figure 9')]:
        add(gene,'Karaoglu2026',locator+'; Methods: in vivo transplantation',
            modification=f'CRISPR-Cas9 {gene} guide; stable selected population',editing_verification='T7E1 and TIDE; protein-level depletion not demonstrated',
            species=species,cell_model=model,challenge_class=challenge,effector_source='immunocompetent CD1 host',matched_control='paired opposite-flank non-targeting graft',
            time_relative_to_challenge='day 0 through 14' if challenge=='allogeneic' else 'day 0 through 7',
            biological_n='4 mice per genetic condition',technical_replication='paired grafts per mouse',immune_validation='separately engineered spheroid grafts',
            protection_result='no durable survival advantage; all grafts cleared',uncertainty='no survival equivalence interval; figure caption does not declare uncertainty convention',
            limitations='Post-immune secretion not measured. Allogeneic text/figure differs on exact clearance day (7 versus 14); xenograft cell dose conflicts (methods 5M versus caption 0.5M).')

for gene,line,control in [('CXCL10','C10G','WT HUES8'),('STAT1','ST1L','GAPDH-luciferase HUES8')]:
    add(gene,'Sintov2022','Supplemental Figure S4F; mmc1.pdf pages 8-9,21',
        modification=f'clonal HDR knockout {line}',editing_verification='genotyping/sequencing; protein and IFNg response checks (Figure 4)',
        species='human',cell_model='SC-islet grafts in NSG mice',challenge_class='none',effector_source='none reported for this baseline assay',
        matched_control=control,time_relative_to_challenge='12-15 weeks after transplantation; no PBMC exposure stated for S4F',
        biological_n='3 plotted observations per group; caption does not explicitly give biological n for panel F',technical_replication='not specified',
        basal_glucose='overnight fasting (16h)',stimulated_glucose='2 g/kg injected glucose; blood at 30 minutes',
        basal_insulin='ratio denominator; absolute concentrations not supplied in S4F',stimulated_insulin='ratio numerator; absolute concentrations not supplied in S4F',
        absolute_stimulated_output='unavailable in S4F',stimulation_index='post-glucose / fasting human insulin; graph only',
        normalization_denominator='same-assay fasting human insulin',uncertainty='mean +/- SD; no retention margin/interval',
        control_responsiveness='ratios above 1 visually; small group and no formal responsiveness test',function_state='F1',
        adverse_function_flag='not established',immune_validation='none in this assay',protection_result='not applicable',
        limitations='Comparator only, never added to primary cohort. Baseline graft maturation assay cannot be re-labelled post-PBMC GSIS.')
    add(gene,'Sintov2022','Figure 4G-J; Supplemental Figure S4G-I',
        modification=f'clonal HDR knockout {line}',editing_verification='Figure 4; Supplemental Figure S4A',species='human',
        cell_model='SC-islets in allogeneic immune-cell co-culture',challenge_class='allogeneic',effector_source='human PBMCs; purified T/NK cells in distinct subexperiments',
        matched_control=control,time_relative_to_challenge='48 hour co-culture',
        biological_n='6 PBMC donors; 2 T-cell donors; 4 NK donors in main Figure 4G (supplement says 2 NK donors)' if gene=='CXCL10' else '2 PBMC/NK donors per Figure 4H',
        technical_replication='n=4 per PBMC/NK donor; n=2-3 per T-cell donor (source notation)' if gene=='CXCL10' else 'n=4 per donor (source notation)',
        immune_validation='separately engineered clone',protection_result='improved viability against PBMC/T/NK challenges' if gene=='CXCL10' else 'pooled-screen rescue not reproduced; more apoptosis after NK co-culture',
        uncertainty='mean +/- SD; source group t tests',limitations='Apoptosis is not secretion. STAT1 is a pooled-versus-clonal counterexample.')
add('CXCL10','Sintov2022','Figure 5A-B; Supplemental Figure S5A; Methods',
    modification='C10G clonal CXCL10 knockout',editing_verification='genotyping; CXCL10 protein/secretion depletion',species='human',
    cell_model='SC-islet kidney grafts in MHC-null NSG mice',challenge_class='allogeneic',effector_source='PBMCs from two human donors',
    matched_control='WT grafts with PBMCs; each genotype without PBMCs',time_relative_to_challenge='week 11 and serial measurements through week 17 after PBMCs',
    biological_n='6-7 mice/group with PBMC; 3/group without; 2 PBMC donors',technical_replication='not specified',
    basal_glucose='overnight fasting before glucose injection',basal_insulin='contemporaneous fasting values not presented in plotted time course',
    stimulated_glucose='2 g/kg glucose; human insulin sampled at 30 minutes',stimulated_insulin='fold relative to pre-PBMC post-glucose measurement',
    absolute_stimulated_output='not reported as absolute concentration in Figure 5B/S5A',stimulation_index='not calculated from contemporaneous basal and stimulated samples in reported time course',
    normalization_denominator='pre-PBMC post-glucose insulin (time 0)',uncertainty='mean +/- SD; no noninferiority margin or retention interval',
    control_responsiveness='post-glucose output present in no-PBMC controls; glucose response itself not established by this plot',
    function_state='F2',adverse_function_flag='none established',immune_validation='separate clone transplanted and challenged',
    protection_result='greater persistence and post-glucose endocrine output than challenged WT',
    limitations='F2 rather than F3/F4: post-glucose timecourse is not paired GSIS. Source text reports endpoint week 17; Figure 5C-D says week 18.')

original_fields=fields.copy()
assert len(rows)==32
assert sum(r['primary_candidate']=='true' for r in rows)==27
assert sum(r['primary_candidate']=='false' for r in rows)==5
assert sum(r['function_state']=='F0' and r['primary_candidate']=='true' for r in rows)==20
assert sum(r['function_state']=='F1' and r['primary_candidate']=='true' for r in rows)==7

# Source-resolved baseline observations are distinct from post-immune retention.
# No numerical value is digitized from a plot. FE018 clarifications are documented
# in evidence/functional_revision3_corrections.json against the preserved v2 row.
baseline_specs={
 'FE017':('Mouse primary islets','Cai 2020','ED2c','Mean increase; limited','Reduced after islet dispersion','Concentration / DNA','No significant edit-NT difference','Baseline; 72 h after transduction','data/raw/literature/uniform20261025/Cai2020_ED2.jpg'),
 'FE018':('Human SC-beta cells','Cai 2020','ED10c','Increase in WT and KO','Responsive','Concentration / 1,000 beta cells','Source reports no impairment','Baseline; no immune challenge','data/raw/literature/uniform20261025/Cai2020_ED10.jpg'),
 'FE019':('Human planar SC-beta cells','MacDonald 2025','ED7D','Absent in WT and KO','No GSIS; KCl responsive','Amount / 100,000 cells','Neither group glucose-responsive','Baseline; no immune challenge','data/raw/literature/rnls2025/mmc1_page7.png'),
 'FE020':('Mouse beta-TC-6 spheroids','Karaoglu 2026','Fig6A','Little/no glucose response','Ratio near one','Ratio only','No significant edit-NT difference','Baseline; pretransplant','data/raw/literature/karaoglu2026/fimmu-17-1759835-g006.jpg'),
 'FE021':('Human EndoC-betaH1 spheroids','Karaoglu 2026','Fig6B','Response observed','Responsive','Ratio only','No significant edit-NT difference','Baseline; pretransplant','data/raw/literature/karaoglu2026/fimmu-17-1759835-g006.jpg'),
 'FE024':('Mouse beta-TC-6 spheroids','Karaoglu 2026','Fig6A','Edited ratio below one','Ratio near one','Ratio only','Lower ratio vs NT; weak control','Baseline; pretransplant','data/raw/literature/karaoglu2026/fimmu-17-1759835-g006.jpg'),
 'FE025':('Human EndoC-betaH1 spheroids','Karaoglu 2026','Fig6B','Response observed','Responsive','Ratio only','No significant edit-NT difference','Baseline; pretransplant','data/raw/literature/karaoglu2026/fimmu-17-1759835-g006.jpg'),
}
baseline=[]
for row in rows:
    if row['record_id'] not in baseline_specs: continue
    model,study,locator,response,control,output,comparison,timing,figure=baseline_specs[row['record_id']]
    baseline.append(dict(record_id=row['record_id'],gene=row['gene'],model_short=model,
        study_short=study,source_locator_short=locator,observed_glucose_response=response,
        control_response=control,output_reporting=output,comparison_observation=comparison,
        immune_challenge_timing=timing,postimmune_function_assessed='no',
        source_figure_file=figure,source_file=row['source_file'],source_url=row['source_url'],
        biological_n=row['biological_n'],technical_replication=row['technical_replication'],
        observation_basis='Source plot/caption and Results; qualitative observation only, no bar digitization',
        comparison_limit='A baseline observation is credited as observed; formal equivalence/noninferiority is not established and no post-immune secretion was measured',
        source_limitations=row['limitations']))
    baseline[-1].update(figure_identity=row['gene']+' | '+model+'\n'+study+' '+locator,
        figure_response=response+'\nControl: '+control,
        figure_output=output,
        figure_comparison=comparison+'\nNo immune challenge')
assert len(baseline)==7
with (ROOT/'evidence'/'baseline_experiments.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(baseline[0]),delimiter='\t');w.writeheader();w.writerows(baseline)

# Append literature-coverage records, not additional experiments or scored candidates.
additional={'LGSN':'Lgsn','PLIN4':'Plin4','TGM6':'Tgm6','GRM2':'Grm2','Or8b47':'Olfr909'}
discovery={**mouse,**additional}
for gene,original in additional.items():
    add(gene,'Cai2020','Figure 1 (discovery identity only); revision search 2026-10-25',
        modification=f'coverage of original discovery label {original}; no new experiment',
        editing_verification='not applicable to a literature-coverage record',
        species='mouse discovery; human and rodent literature searched',
        cell_model='NIT-1 discovery context; no qualifying exact-edit secretion assay located',
        challenge_class='literature coverage: baseline, autoimmune, allogeneic',
        time_relative_to_challenge='no qualifying direct secretion comparison located',
        function_state='not applicable: coverage record',
        immune_validation='original pooled discovery selection only',
        protection_result='no additional independent beta-cell immune-protection validation located',
        limitations='Not an experimental comparison. Search absence is bounded, not proof of biological absence. Quantitative cohort unchanged.' +
                    (' Or8b47 is the current mouse symbol (NCBI Gene 258873); no human counterpart is assigned.' if gene=='Or8b47' else ''))

# A completed source assessment does not add another extracted experiment.
add('RNLS','Bode2023','Full Methods/Results/Discussion; main Figures 1-7 legends; official Supplemental Figures 1-5 and legends; PMID37216639',
    source_file='data/raw/literature/bode2023_final/PMC10382656_article_AX.txt',
    modification='source assessment of Rnls-mutant NIT-1 cells; not a single extracted comparison',
    editing_verification='study methods report deep-sequencing and quantitative-PCR verification of mutant line',
    species='mouse',cell_model='NIT-1 cells and paired subcutaneous grafts in NOD.scid recipients',
    challenge_class='source covers autoimmune exposure and unchallenged controls',
    effector_source='autoreactive splenocytes from diabetic NOD mice; separate T-cell coculture',
    matched_control='study uses WT/non-targeting cells and paired grafts; no single contrast assigned to source-assessment record',
    time_relative_to_challenge='multiple source assays: baseline flux, 72-hour T-cell coculture, and graft retrieval after splenocyte transfer; no secretion timing',
    biological_n='not applicable to aggregate source-assessment record; experiment-specific replication remains in source legends',
    technical_replication='not applicable to aggregate source-assessment record',
    function_state='not applicable: completed full-text source assessment',
    immune_validation='source reports paired graft weights, infiltrating immune-cell phenotyping and local immune modulation',
    protection_result='source-level context: Rnls-mutant graft weight advantage under autoimmune exposure (Figures 2B-C and 4A); not a secretion endpoint',
    limitations='Full public primary text and all six official supplemental PDFs assessed. No direct insulin/C-peptide secretion, GSIS or glucose-function endpoint located. Basal ECAR/OCR, graft weight and immune-cell function do not establish beta-cell secretion. This source-assessment row does not enumerate or pool the paper\'s immune experiments.')
add('RNLS','Bode2024','Figure 1B-I; Methods: beta-cell transplantation',
    modification='Rnls-mutant NIT-1 cells; original mutant line',editing_verification='line described previously; source cites Cai2020',
    species='mouse',cell_model='paired WT/Rnls-mutant subcutaneous NIT-1 grafts in NOD.scid or NSG recipients',
    challenge_class='autoimmune',effector_source='autoreactive NOD splenocytes with or without NKp46-positive cell depletion',
    matched_control='paired WT graft and undepleted splenocyte condition',
    time_relative_to_challenge='grafts retrieved17 days after splenocyte transfer; transfer4 days after grafting',
    biological_n='9 control paired graft mice and 5 NKp46-depleted paired graft mice across stated experiments',
    technical_replication='paired grafts; source Figure1C control group includes4 prior2023 mice and5 new mice',
    immune_validation='graft weight and immune-cell phenotyping; conditional NKp46 depletion test',
    protection_result='Rnls-mutant graft advantage depends on NKp46-positive innate lymphoid cells in this model',
    uncertainty='source mean and paired two-tailed tests; no secretion endpoint interval',
    limitations='No direct insulin-secretion assay in located main text/supplement. Four of nine control mice and scRNAseq data overlap Bode2023; not wholly independent replication.')
add('RNLS','Bode2024','Figure 4H-I; Methods: SCBC/PBMC co-culture; Supplemental Figure S4',
    modification='RNLS-mutant human stem-cell-derived beta-like cells',editing_verification='mutant line from prior study; differentiation identity via C-peptide/NKX6.1 staining',
    species='human',cell_model='differentiated SC-beta-like cell culture; survival gate is alive CD45-negative cells, not pure beta-marker-gated population',
    challenge_class='allogeneic in vitro',effector_source='HLA-mismatched human PBMCs, IL2 stimulated',
    matched_control='WT differentiated culture in matched co-culture; without-PBMC conditions shown',
    time_relative_to_challenge='72-hour co-culture',
    biological_n='one of two independent experiments shown; n=2-3 per source caption; donor/batch allocation not fully resolved',
    technical_replication='n=2-3 reported; not assumed to be independent immune donors',
    immune_validation='separately engineered cells in short HLA-mismatched co-culture',
    protection_result='favorable alive CD45-negative cell count and alive/dead ratio versus challenged WT',
    uncertainty='mean +/- SD; source Figure4H one-way ANOVA/Dunnett and Figure4I two-tailed t-test; no retention interval',
    limitations='Favorable short in-vitro survival observation, not durable graft or clinical protection. C-peptide/NKX6.1 stains are identity, not secreted-insulin measurements; no GSIS endpoint in located main text/supplement.')

fields += ['record_role','discovery_member','original_mouse_label','quantitative_eligible',
           'experiment_group','measurement_timing','basal_stimulated_assessment',
           'absolute_output_assessment','control_response_assessment','comparative_function_result',
           'function_measurement_status','immune_protection_status','search_date',
           'observed_glucose_response','baseline_comparison_observation','immune_challenge_during_function_assay','baseline_figure_file','source_access_status']
for row in rows:
    number=int(row['record_id'][2:])
    coverage=33<=number<=37
    source_assessment=number==38
    if source_assessment:
        for key,value in row.items():
            if value==N:row[key]='not applicable: no direct secretion assay located in assessed full text and supplement'
        row['retention_margin']='not applicable: no secretion-retention analysis located in assessed source'
        row['retention_estimate_interval']='not applicable: no secretion-retention estimate located in assessed source'
        row['adverse_function_flag']='not assessed: absence of a secretion endpoint is not adverse function'
    row.update(record_role='supplemental_literature_coverage' if coverage else
               'completed_full_text_source_assessment' if source_assessment else
               'revision_primary_candidate_comparison' if number>=39 else
               'primary_candidate_comparison' if row['primary_candidate']=='true' else 'interpretation_control_comparison',
               discovery_member=str(row['gene'] in discovery).lower(),
               original_mouse_label=discovery.get(row['gene'],'not in original discovery set'),
               quantitative_eligible=row['primary_candidate'],
               search_date='2026-10-25' if coverage else '2026-10-24')
    row['experiment_group']=('not an extracted experiment' if source_assessment else 'not an experiment' if coverage else
                             'Cai2020_pooled_discovery' if int(row['record_id'][2:])<=12 and row['study']=='Cai2020' else
                             'Sintov2022_pooled_allogeneic_screen' if int(row['record_id'][2:])<=12 else
                             row['study']+'_'+row['source_locator'].split(';')[0])
    row['measurement_timing']='not located' if coverage else row['time_relative_to_challenge']
    row['basal_stimulated_assessment']=('not located' if coverage else
        'paired basal and stimulated assay reported' if row['function_state']=='F1' else
        'post-glucose output; contemporaneous basal not presented' if row['function_state']=='F2' else
        'candidate-specific direct secretion not measured/reported in this comparison')
    row['absolute_output_assessment']=row['absolute_stimulated_output']
    row['control_response_assessment']=row['control_responsiveness']
    row['comparative_function_result']=('no qualifying comparison located' if coverage else row['adverse_function_flag'])
    row['function_measurement_status']=('not located' if coverage else
        'measured but inconclusive for retained function' if row['function_state'] in ['F1','F2'] else
        'not measured/reported in located comparison')
    row['immune_protection_status']=row['protection_result']
    details=next((b for b in baseline if b['record_id']==row['record_id']),None)
    row['observed_glucose_response']=details['observed_glucose_response'] if details else ('no glucose-response secretion assay located in assessed source' if source_assessment else 'not measured in this comparison')
    row['baseline_comparison_observation']=details['comparison_observation'] if details else 'not applicable'
    row['immune_challenge_during_function_assay']='none; baseline' if details else ('post-immune output assay' if row['function_state']=='F2' else 'no direct function assay extracted')
    row['baseline_figure_file']=details['source_figure_file'] if details else 'not applicable'
    row['source_access_status']='complete public primary text and all six official supplemental PDFs assessed' if source_assessment else 'primary full text/supplement assessed' if not coverage else 'bounded literature-coverage record'
    if number>=38:row['search_date']='2026-10-25'
    if source_assessment:
        row['measurement_timing']='not applicable: source-level assessment; no direct secretion timing'
        row['basal_stimulated_assessment']='no direct basal/stimulated secretion comparison located in assessed full text and supplement'
        row['function_measurement_status']='no direct secretion endpoint located in assessed full text and supplement'

for suffix,delimiter in [('csv',','),('tsv','\t')]:
    p=ROOT/'evidence'/f'functional_evidence.{suffix}'
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter=delimiter);w.writeheader();w.writerows(rows)
shutil.copyfile(ROOT/'evidence'/'functional_evidence.tsv', ROOT/'functional_evidence.tsv')
assert all(any(r['gene']==g for r in rows) for g in candidates)
assert not any(r['function_state'] in ['F3','F4'] for r in rows)
assert all(r['joint_gate']=='not established' for r in rows)
summaries=[]
matrix=[]
for gene in discovery:
    baseline_ids=[r['record_id'] for r in rows if r['gene']==gene and r['function_state']=='F1']
    locators='; '.join(dict.fromkeys(r['study']+' '+r['source_locator'].split(';')[0] for r in rows if r['record_id'] in baseline_ids))
    auto_locators='Cai2020 Fig1' + ('; Fig2a-l' if gene=='RNLS' else '')
    allo_locators=('Sintov2022 Fig1G-H (pooled)' if gene in candidates else 'revision20261025 search')
    if gene=='RNLS': allo_locators+='; Cai2020 ED3b-c; Karaoglu2026 Fig8; Bode2024 Fig4H-I'
    if gene=='HIVEP2': allo_locators+='; Karaoglu2026 Fig8'
    matrix.append(dict(gene=gene,original_mouse_label=discovery[gene],quantitative_eligible=str(gene in candidates).lower(),
        analysis_group='six scored counterparts' if gene in candidates else 'supplemental discovery literature coverage',
        baseline_state='M' if baseline_ids else 'NL',
        baseline_text='Measured; see model-study panel' if baseline_ids else 'No qualifying assay located within specified search',
        baseline_timing='unchallenged / pretransplant' if baseline_ids else 'not located',
        baseline_basal_stimulated='reported' if baseline_ids else 'not located',
        baseline_absolute_output='ratios only' if gene=='HIVEP2' else 'some normalized absolute plots; numeric extraction absent' if gene=='RNLS' else 'not located',
        baseline_control_responsiveness='human responsive; mouse weak/unresponsive; see model-study rows' if gene=='HIVEP2' else 'Cai human SC-beta responsive; MacDonald planar SC-beta unresponsive; other model-specific rows separate' if gene=='RNLS' else 'not located',
        baseline_comparative_result='observed baseline responses and edited-control comparisons separated in baseline_experiments.tsv' if baseline_ids else 'no qualifying comparison located',
        baseline_model_modifier='mouse beta-TC-6 and human EndoC-betaH1 spheroids' if gene=='HIVEP2' else 'mouse primary islets/spheroids; human SC-beta cells/spheroids' if gene=='RNLS' else 'exact loss-of-function beta-cell secretion required',
        baseline_locators=locators if baseline_ids else 'evidence/literature_search.json',
        baseline_record_ids=','.join(baseline_ids),
        post_autoimmune_state='NM',post_autoimmune_text='Not measured for individual edit',
        post_autoimmune_timing='after immune exposure',post_autoimmune_basal_stimulated='not measured/reported',
        post_autoimmune_absolute_output='not measured/reported',post_autoimmune_control_responsiveness='not assessed',
        post_autoimmune_comparative_result='not estimable',post_autoimmune_model_modifier='mouse subcutaneous NIT-1 pooled discovery; separate NIT-1 and kidney-capsule primary-islet RNLS validation' if gene=='RNLS' else 'mouse subcutaneous NIT-1 pooled discovery',
        post_autoimmune_locators=auto_locators,
        post_allogeneic_state='NM' if gene in candidates else 'NL',
        post_allogeneic_text='Not measured for individual edit' if gene in candidates else 'No qualifying assay located',
        post_allogeneic_timing='after immune exposure' if gene in candidates else 'not located',
        post_allogeneic_basal_stimulated='not measured/reported' if gene in candidates else 'not located',
        post_allogeneic_absolute_output='bulk pooled-graft insulin is not edit-specific' if gene in candidates else 'not located',
        post_allogeneic_control_responsiveness='not established per edit' if gene in candidates else 'not located',
        post_allogeneic_comparative_result='not estimable',
        post_allogeneic_model_modifier='human SC-islet pooled screen; independent mouse grafts; human-to-mouse xenografts kept separate' if gene in ['HIVEP2','RNLS'] else 'human SC-islet pooled screen' if gene in candidates else 'no qualifying independent assay located',
        post_allogeneic_locators=allo_locators,
        immune_protection_state='MIXED' if gene=='RNLS' else 'AC' if gene=='HIVEP2' else 'SCREEN',
        immune_protection_text='Autoimmune validation; favorable short human PBMC co-culture survival; no durable mouse allo protection' if gene=='RNLS' else 'No durable protection in tested mouse allograft' if gene=='HIVEP2' else 'Pooled discovery only; no independent validation located within specified search',
        immune_protection_locators='Cai2020 Fig2,ED3; Bode2024 Fig1,4H-I; Karaoglu2026 Fig8' if gene=='RNLS' else 'Karaoglu2026 Fig8' if gene=='HIVEP2' else 'Cai2020 Fig1',
        figure_protection_text='Auto validated: C Fig2; B Fig1\n72 h PBMC benefit: B Fig4\nNo durable mouse allo: C ED3; K Fig8' if gene=='RNLS' else 'No durable mouse allo\nK Fig8' if gene=='HIVEP2' else 'Pooled nomination\nC Fig1',
        adequate_postimmune_function='not established',joint_result='not established',
        interpretation='Evidence coverage; not a biological failure rate or gene-wide function grade'))
    summaries.append(dict(gene=gene,original_mouse_label=discovery[gene],quantitative_eligible=str(gene in candidates).lower(),
        baseline_secretion_measurement='direct baseline assay; observed response and comparison reported separately by model-study' if baseline_ids else 'no qualifying direct assay located within specified search',
        post_autoimmune_secretion='not measured/reported for individual edit in located immune comparison',
        post_allogeneic_secretion='not measured/reported for individual edit in located immune comparison' if gene in candidates else 'no qualifying assay located',
        independent_autoimmune_validation='supported in Cai mouse models' if gene=='RNLS' else 'not established; discovery-screen selection only',
        independent_allogeneic_validation='no durable protection in tested mouse spheroid model' if gene=='HIVEP2' else 'favorable72-hour human PBMC co-culture survival; no durable protection in tested mouse graft models' if gene=='RNLS' else 'not located within specified search',
        both_challenges_with_retained_function='not established',
        caveat='Evidence coverage, not a biological failure rate. No immutable gene grade. Five supplemental targets do not enter quantitative test.'))
with (ROOT/'evidence'/'functional_candidate_summary.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter='\t');w.writeheader();w.writerows(summaries)
with (ROOT/'evidence'/'functional_matrix.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(matrix[0]),delimiter='\t');w.writeheader();w.writerows(matrix)
with (ROOT/'evidence'/'functional_matrix.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(matrix[0]));w.writeheader();w.writerows(matrix)
controls=[dict(gene='CXCL10',scope='preclinical interpretation control; outside original eleven',baseline='ratio GSIS before PBMC challenge; retention inconclusive (S4F)',post_allogeneic='post-glucose output relative to pre-PBMC; basal/stimulated pairing absent (Fig5B/S5A)',protection='separate clone protection in allogeneic coculture/graft',clinical_endpoint='not measured',retention='not established'),
          dict(gene='STAT1',scope='preclinical interpretation control; outside original eleven',baseline='ratio GSIS before PBMC challenge; retention inconclusive (S4F)',post_allogeneic='secretion not measured in clonal immune coculture (Fig4H)',protection='pooled hit not reproduced; greater NK apoptosis',clinical_endpoint='not measured',retention='not established')]
with (ROOT/'evidence'/'functional_interpretation_controls.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(controls[0]),delimiter='\t');w.writeheader();w.writerows(controls)
counts=dict(original_comparisons=32,primary_candidate_comparisons=27,interpretation_control_comparisons=5,
            revision_primary_candidate_comparisons=2,completed_full_text_source_assessment_records=1,
            total_extracted_comparisons=34,original_studies=4,new_primary_followup_studies=2,baseline_model_study_comparisons=7,
            supplemental_literature_coverage_records=5,total_records=len(rows),original_primary_F0=20,original_primary_F1=7,
            original_primary_F2_to_F4=0,scored_candidates=6,discovery_targets_covered=11,
            baseline_direct_measurement_scored_candidates=2,baseline_direct_measurement_all_discoveries=2,
            adequate_postimmune_assessment_scored_candidates=0,adequate_postimmune_assessment_all_discoveries=0,
            zero_of_six_interpretation='evidence coverage; no tested denominator for biological preserved-function rate')
(ROOT/'evidence'/'functional_coverage_counts.json').write_text(json.dumps(counts,indent=2)+'\n',encoding='utf-8')
print('Wrote 40 records: 34 extracted comparisons, 5 coverage records, 1 completed full-text source assessment; 11-target matrix and 7-row baseline panel; quantitative cohort fixed at six.')
