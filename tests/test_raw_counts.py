import unittest
import numpy as np
import pandas as pd
from src.raw_counts import cpm, contrast, theta

class RawCountUnitTests(unittest.TestCase):
    def test_normalization_is_invariant_to_sequencing_depth(self):
        a=pd.DataFrame({'mouse1':[3.,7.],'mouse2':[6.,14.]})
        n=cpm(a)
        np.testing.assert_allclose(n.mouse1,n.mouse2)
        np.testing.assert_allclose(n.sum(),[1e6,1e6])

    def test_zero_library_is_rejected(self):
        with self.assertRaises(ValueError): cpm(pd.DataFrame({'mouse1':[0.,0.]}))

    def test_mouse_means_not_guide_pseudoreplicates(self):
        # One row per gene after guide averaging, with three real sample columns.
        g=pd.DataFrame({'treated1':[1.],'treated2':[5.],'control1':[2.]})
        self.assertEqual(contrast(g,['treated1','treated2'],['control1']).iloc[0],1.)
        self.assertEqual(contrast(g,['treated2'],['control1']).iloc[0],3.)

    def test_theta_ties_against_direct_pairs(self):
        v=np.array([0.,1.,1.,2.]);selected=np.array([False,True,True,False])
        a=v[selected];b=v[~selected]
        direct=np.mean((a[:,None]>b[None,:])+.5*(a[:,None]==b[None,:]))
        self.assertEqual(theta(v,selected),direct)

    def test_direction_and_background_exclusion(self):
        self.assertEqual(theta(np.array([1.,2.,3.]),np.array([False,False,True])),1.)
        self.assertEqual(theta(np.array([1.,2.,3.]),np.array([True,False,False])),0.)

if __name__=='__main__': unittest.main()
