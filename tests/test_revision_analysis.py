"""Independent small-universe checks and revision statistical contracts."""
from itertools import combinations,product
from pathlib import Path
import json
import math
import unittest
import numpy as np
import pandas as pd
from src.revision_analysis import (theta_from_rank_sums,uniform_theta_draws,
    no_hit_probability,stratified_theta,matched_randomization,weighted_rank_tuples)
from src.transfer import rankdata,superiority

ROOT=Path(__file__).resolve().parents[1]


class RevisionAnalysisTests(unittest.TestCase):
    def test_rank_sum_matches_brute_force_with_ties(self):
        scores=np.array([0.,1.,1.,3.,4.])
        for indices in combinations(range(len(scores)),2):
            a=scores[list(indices)];b=np.delete(scores,indices)
            direct=np.mean((a[:,None]>b[None,:])+.5*(a[:,None]==b[None,:]))
            self.assertAlmostEqual(theta_from_rank_sums(rankdata(scores)[list(indices)].sum(),5,2),direct)

    def test_uniform_small_universe_exact_mean_variance(self):
        scores=np.arange(7,dtype=float)
        exact=np.array([superiority(scores,np.array(c)) for c in combinations(range(7),3)])
        simulated=uniform_theta_draws(scores,3,80_000,321)
        self.assertAlmostEqual(exact.mean(),.5)
        self.assertLess(abs(simulated.mean()-exact.mean()),.003)
        self.assertLess(abs(simulated.var()-exact.var()),.002)

    def test_stratified_statistic_and_sampler_against_exact(self):
        frame=pd.DataFrame({'gene':['a','b','c','d','e','f'],
            'published_score':[0.,1.,3.,0.,2.,2.],'guide_count':[4]*6,
            'control_quintile':[0]*3+[1]*3,'pretransplant_quintile':[0]*6}).set_index('gene')
        selected=np.array([True,False,False,False,True,False])
        observed,draws,records=matched_randomization(frame,['a','e'],80_000,444)
        expected=stratified_theta(frame.published_score.to_numpy(),['x']*3+['y']*3,selected)
        self.assertAlmostEqual(observed,expected)
        exact=[]
        for i,j in product(range(3),range(3,6)):
            sel=np.zeros(6,dtype=bool);sel[[i,j]]=True
            exact.append(stratified_theta(frame.published_score.to_numpy(),['x']*3+['y']*3,sel))
        self.assertLess(abs(draws.mean()-np.mean(exact)),.003)
        self.assertLess(abs(draws.var()-np.var(exact)),.002)
        self.assertTrue(records.candidate_count.eq(1).all())

    def test_weighted_rank_rejection_preserves_declared_conditional_law(self):
        n=5;alphas=np.array([2.,1.]);rng=np.random.Generator(np.random.PCG64(910))
        sampled,rejected=weighted_rank_tuples(rng,n,alphas,100_000)
        self.assertTrue((sampled[:,0]!=sampled[:,1]).all());self.assertGreater(rejected,0)
        weights=np.stack([np.diff((np.arange(n+1)/n)**a) for a in alphas])
        law=weights[0,:,None]*weights[1,None,:];np.fill_diagonal(law,0);law/=law.sum()
        actual=np.bincount(sampled[:,0]*n+sampled[:,1],minlength=n*n).reshape(n,n)/len(sampled)
        self.assertLess(np.abs(actual-law).max(),.004)

    def test_hypergeometric_no_hit_identity_and_boundaries(self):
        self.assertAlmostEqual(no_hit_probability(10,3,2),math.comb(7,2)/math.comb(10,2))
        self.assertEqual(no_hit_probability(10,0,6),1.)
        self.assertEqual(no_hit_probability(10,10,1),0.)
        self.assertEqual(no_hit_probability(10,10,0),1.)

    def test_shared_constant_preserves_all_orderings_and_theta(self):
        x=np.array([-2.,-.2,0.,0.,.7,2.]);ids=np.array([1,3])
        shifted=x+1.234
        np.testing.assert_equal(rankdata(x),rankdata(shifted))
        self.assertEqual(superiority(x,ids),superiority(shifted,ids))
        self.assertNotEqual(int((x>0).sum()),int((shifted>0).sum()))

    def test_calibration_contract_and_original_result_unchanged(self):
        r=json.loads((ROOT/'results/revision_summary.json').read_text())
        primary=pd.read_csv(ROOT/'results/transfer_summary.tsv',sep='\t').iloc[0]
        c=r['primary_null_calibration']
        self.assertAlmostEqual(c['observed_theta'],primary.theta,places=14)
        self.assertEqual(c['exceedances'],int(primary.null_exceedances))
        self.assertAlmostEqual(c['p_upper_tail_plus_one'],primary.p_one_sided,places=14)
        self.assertLess(c['central_95_lower'],c['observed_theta'])
        self.assertGreater(c['central_95_upper'],c['observed_theta'])
        self.assertLessEqual(c['reference_null_rejection_fraction'],.05)
        h=r['no_hit_calibration'];self.assertEqual(h['positive_genomewide_supported_genes'],98)
        self.assertAlmostEqual(h['uniform_random_set_probability_zero_hits'],.9507286242496122,places=14)
        matched=r['secondary_joint_matched_test']
        self.assertEqual(matched['common_universe_n'],11474)
        self.assertEqual(matched['universe_missing_covariates'],215)
        self.assertEqual(matched['n_occupied_strata'],6)
        individual=pd.read_csv(ROOT/'results/raw_matched_background.csv')
        self.assertAlmostEqual(matched['observed_T'],individual.candidate_vs_matched_probability.mean(),places=14)
        self.assertTrue(r['common_universe_constant_shift']['rank_vectors_identical'])
        self.assertTrue(r['common_universe_constant_shift']['theta_identical'])


if __name__=='__main__': unittest.main()
