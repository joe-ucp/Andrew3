"""Scientific contract and frozen-provenance checks over generated artifacts."""
from pathlib import Path
import csv
import hashlib
import json
import unittest
import numpy as np
import pandas as pd
from src.source_inventory import inventory,relative_path,inventory_path,validate_source_locator

ROOT=Path(__file__).resolve().parents[1]
CANDIDATES={'HIVEP2','RNLS','GCN1','MED11','ZBED3','ZNF28'}

class IntegrationTests(unittest.TestCase):
    def test_every_frozen_source_hash_and_inventory_row(self):
        expected=inventory(ROOT)
        with inventory_path(ROOT).open(encoding='utf-8',newline='') as f:
            actual=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual({r['filename'] for r in actual},{r['filename'] for r in expected})
        e={r['filename']:r for r in expected}
        for row in actual:
            self.assertEqual(row['sha256'],e[row['filename']]['sha256'])
            self.assertTrue(row['retrieval_utc'] and row['url_or_accession'] and row['access_note'] and row['role'])
            self.assertNotIn('_UCP',row['filename'])

    def test_source_paths_are_portable_and_confined(self):
        self.assertEqual(relative_path('data\\raw\\a.txt'),'data/raw/a.txt')
        for value in ['../secret','/tmp/secret','C:/private/file']:
            with self.assertRaises(ValueError): relative_path(value)

    def test_complete_discovery_denominator_and_no_imputed_scores(self):
        m=pd.read_csv(ROOT/'candidate_mapping.tsv',sep='\t')
        self.assertEqual(len(m),11);self.assertTrue(m.discovery_label.is_unique)
        self.assertEqual(set(m.loc[m.included,'human_symbol']),CANDIDATES)
        self.assertEqual(int((~m.included).sum()),5)
        self.assertTrue(m.loc[~m.included,'exclusion_reason'].notna().all())
        d=pd.read_csv(ROOT/'data/processed/allogeneic_scores.csv')
        self.assertEqual(len(d),11689);self.assertTrue(d.gene.is_unique)
        self.assertTrue(np.isfinite(d.published_score).all())
        np.testing.assert_allclose(d.reported_q,10.**(-d.published_neg_log10_fdr),rtol=1e-12,atol=0.)

    def test_primary_regression_and_interpretation_contract(self):
        t=pd.read_csv(ROOT/'results/transfer_summary.tsv',sep='\t')
        p=t.loc[t.analysis.eq('primary')].iloc[0]
        self.assertEqual(p.n_universe,11689);self.assertEqual(p.n_candidates,6)
        self.assertEqual(p.draws,1_000_000);self.assertEqual(p.seed,20261025)
        self.assertAlmostEqual(p.theta,0.42096778795400724,places=14)
        self.assertEqual(p.null_exceedances,743859)
        self.assertAlmostEqual(p.p_one_sided,(743859+1)/(1_000_000+1),places=14)
        self.assertEqual(p.positive_scores,2);self.assertEqual(p.positive_supported,0)
        self.assertTrue((t.theta<.5).all())

    def test_functional_comparisons_do_not_invent_joint_support(self):
        f=pd.read_csv(ROOT/'functional_evidence.tsv',sep='\t',keep_default_na=False)
        self.assertTrue(f.record_id.is_unique)
        flags=f.primary_candidate.astype(str).str.lower().eq('true')
        self.assertEqual(set(f.loc[flags,'gene']),CANDIDATES)
        primary=f.record_role.eq('primary_candidate_comparison')
        controls=f.record_role.eq('interpretation_control_comparison')
        coverage=f.record_role.eq('supplemental_literature_coverage')
        revision=f.record_role.eq('revision_primary_candidate_comparison')
        assessment=f.record_role.eq('completed_full_text_source_assessment')
        self.assertEqual(len(f),40)
        self.assertEqual(int(primary.sum()),27)
        self.assertEqual(int(controls.sum()),5)
        self.assertEqual(int(coverage.sum()),5)
        self.assertEqual(int(revision.sum()),2)
        self.assertEqual(int(assessment.sum()),1)
        self.assertTrue((primary|controls|coverage|revision|assessment).all())
        self.assertEqual(int((primary|controls|revision).sum()),34)
        self.assertEqual(set(f.loc[controls,'gene']),{'CXCL10','STAT1'})
        supplemental={'LGSN','PLIN4','TGM6','GRM2','Or8b47'}
        self.assertEqual(set(f.loc[coverage,'gene']),supplemental)
        self.assertFalse(flags[coverage].any())
        self.assertTrue(f.loc[coverage,'function_state'].eq('not applicable: coverage record').all())
        self.assertEqual(set(f.loc[primary|controls,'record_id']),{f'FE{i:03d}' for i in range(1,33)})
        self.assertEqual(set(f.loc[revision,'gene']),{'RNLS'})
        self.assertTrue(f.loc[assessment,'function_state'].eq('not applicable: completed full-text source assessment').all())
        self.assertTrue(f.loc[assessment,'experiment_group'].eq('not an extracted experiment').all())
        self.assertTrue(f.loc[assessment,'basal_insulin'].str.contains('no direct secretion assay located').all())
        self.assertTrue(f.loc[assessment,'retention_margin'].str.contains('no secretion-retention analysis located').all())
        self.assertTrue(f.loc[assessment,'source_access_status'].str.contains('all six official supplemental PDFs assessed').all())
        self.assertFalse(f.loc[assessment].map(lambda v: 'unresolved' in str(v) or 'abstract-limited' in str(v)).any().any())
        self.assertEqual(f.loc[primary,'function_state'].value_counts().to_dict(),{'F0':20,'F1':7})
        m=pd.read_csv(ROOT/'evidence/functional_matrix.tsv',sep='\t',keep_default_na=False)
        self.assertEqual(len(m),11);self.assertTrue(m.gene.is_unique)
        self.assertEqual(set(m.gene),CANDIDATES|supplemental)
        self.assertEqual(set(m.loc[m.quantitative_eligible.astype(str).str.lower().eq('true'),'gene']),CANDIDATES)
        self.assertTrue(m.adequate_postimmune_function.eq('not established').all())
        self.assertEqual(int(m.baseline_state.eq('M').sum()),2)
        self.assertEqual(m.loc[m.gene.eq('RNLS'),'immune_protection_state'].item(),'MIXED')
        b=pd.read_csv(ROOT/'evidence/baseline_experiments.tsv',sep='\t',keep_default_na=False)
        self.assertEqual(len(b),7);self.assertTrue(b.record_id.is_unique)
        self.assertTrue(b.postimmune_function_assessed.eq('no').all())
        self.assertIn('Increase in WT and KO',b.loc[b.record_id.eq('FE018'),'observed_glucose_response'].item())
        self.assertIn('Absent in WT and KO',b.loc[b.record_id.eq('FE019'),'observed_glucose_response'].item())
        self.assertIn('1,000 beta cells',b.loc[b.record_id.eq('FE018'),'output_reporting'].item())
        self.assertTrue(f.loc[f.record_id.eq('FE040'),'protection_result'].str.contains('favorable').all())
        self.assertTrue(f.loc[f.record_id.eq('FE040'),'cell_model'].str.contains('not pure beta-marker-gated').all())
        self.assertFalse(f.loc[flags,'function_state'].isin(['F3','F4']).any())
        self.assertTrue(f.joint_gate.eq('not established').all())
        required=['modification','species','cell_model','challenge_class','matched_control',
                  'time_relative_to_challenge','biological_n','technical_replication','basal_glucose',
                  'basal_insulin','stimulated_glucose','stimulated_insulin','absolute_stimulated_output',
                  'normalization_denominator','uncertainty','control_responsiveness','source_locator','source_file']
        self.assertTrue(f[required].map(lambda value:bool(str(value).strip())).all().all())
        referenced_sources=set(f.source_file)|set(b.source_file)|set(b.source_figure_file)
        for source in referenced_sources:
            self.assertGreater(validate_source_locator(ROOT,source)['bytes'],0)

    def test_raw_units_coverage_and_exploratory_separation(self):
        qc=pd.read_csv(ROOT/'results/raw_sample_qc.csv')
        self.assertEqual(sum(qc.sample_class.eq('no_PBMC')),6)
        self.assertEqual(sum(qc.sample_class.eq('hPi')),6)
        r=pd.read_csv(ROOT/'results/raw_candidate_summary.csv')
        self.assertEqual(set(r.gene),CANDIDATES);self.assertTrue(r.guide_count.eq(4).all())
        loo=pd.read_csv(ROOT/'results/raw_leave_one_animal_out.csv')
        self.assertEqual(len(loo),72)
        for _,group in loo.groupby('gene'): self.assertEqual(group.omitted_animal.nunique(),12)
        nt=pd.read_csv(ROOT/'results/raw_exploratory_NT_centered.csv')
        self.assertTrue(nt.NT_guides.eq(941).all())
        self.assertTrue(nt.same_direction_as_published.all())
        self.assertTrue(nt.analysis_status.eq('post_primary_exploratory_different_estimand').all())
        plan=json.loads((ROOT/'protocol/raw_sensitivity.json').read_text())
        self.assertIn('primary_outcome_known',plan['post_primary_outcome_amendment'])

    def test_figures_are_nonempty_and_not_missing(self):
        from PIL import Image
        for name in ['transfer.png','evidence_matrix.png']:
            p=ROOT/'manuscript/figures'/name
            with Image.open(p) as im:
                self.assertGreater(im.width,1000);self.assertGreater(im.height,600)

if __name__=='__main__': unittest.main()
