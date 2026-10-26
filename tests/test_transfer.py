import itertools
import unittest
import numpy as np
from src.transfer import superiority, sample_subsets, competitive_test

class TransferTests(unittest.TestCase):
    def test_ties_against_direct_pairs_all_subsets(self):
        scores=np.array([-2.,0.,0.,3.,7.])
        for m in range(1,5):
            for c in itertools.combinations(range(5),m):
                b=[i for i in range(5) if i not in c]
                direct=np.mean([float(scores[i]>scores[j])+.5*float(scores[i]==scores[j]) for i in c for j in b])
                self.assertAlmostEqual(superiority(scores,c),direct)
    def test_monotone_and_direction(self):
        scores=np.array([-3.,-1.,0.,2.,3.])
        self.assertEqual(superiority(scores,[3,4]),1.)
        self.assertEqual(superiority(scores,[0,1]),0.)
        self.assertAlmostEqual(superiority(scores,[1,3]),superiority(np.exp(scores),[1,3]))
    def test_all_tied(self):
        result=competitive_test(np.ones(8),[1,2],draws=1000)
        self.assertEqual(result['theta'],.5)
        self.assertEqual(result['p_one_sided'],1.)
    def test_sampler_exact_null_small_universe(self):
        rows=sample_subsets(np.random.Generator(np.random.PCG64(6)),4,2,120000)
        self.assertTrue(np.all(rows[:,0]!=rows[:,1]))
        counts={c:0 for c in itertools.combinations(range(4),2)}
        for row in rows: counts[tuple(sorted(row))]+=1
        self.assertTrue(all(abs(n-20000)<600 for n in counts.values()),counts)
        out=competitive_test(np.arange(4.),[2,3],draws=120000,seed=6)
        self.assertAlmostEqual(out['p_one_sided'],1/6,delta=.004)

if __name__=='__main__': unittest.main()
