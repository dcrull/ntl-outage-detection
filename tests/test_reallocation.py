from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from pyproj import Transformer
from rasterio.transform import from_origin
from scipy import sparse

from scripts.reallocation import (upstream_modules, aggregate_complete,
                                  allocate_nightly, native_overlap_matrix)


class ReallocationTests(unittest.TestCase):
    def test_aggregation_requires_full_area_support(self):
        w=sparse.csr_matrix([[.5,.5,0],[0,0,1.]])
        values,support=aggregate_complete(w,np.array([[2,np.nan,0.]]))
        np.testing.assert_equal(values,[np.nan,0.])
        np.testing.assert_allclose(support,[.5,1.])
        values,_=aggregate_complete(w,np.array([[2,4,0.]]))
        np.testing.assert_allclose(values,[3.,0.])

    def test_native_area_weights_recover_constant_and_fail_on_cropped_cell(self):
        cells=pd.DataFrame({'source_row':[945],'source_col':[144]})
        tr=Transformer.from_crs('EPSG:4326','EPSG:32616',always_xy=True)
        x,y=tr.transform(-90+144/240,40-945/240)
        transform=from_origin(x-100,y+100,10,10)
        w,area=native_overlap_matrix(cells,transform,(100,100))
        self.assertGreater(area[0],100000)
        self.assertAlmostEqual(float(w.sum()),1,places=8)
        values,support=aggregate_complete(w,np.full((100,100),7.))
        np.testing.assert_allclose(values,[7.]);np.testing.assert_allclose(support,[1.])
        with self.assertRaisesRegex(ValueError,'does not cover'):
            native_overlap_matrix(cells,transform,(3,3))

    @unittest.skipUnless((Path.home()/'research/ntl-psf-disaggregation/src').exists(),
                         'Original checkout required for reference-operator comparison')
    def test_optimized_allocation_matches_original_with_clouds_and_zeros(self):
        m=upstream_modules(Path.home()/'research/ntl-psf-disaggregation')
        op=m['operator']; kernel=m['kernels'].circular_mean_kernel(radius_m=30,resolution_m=10)
        rng=np.random.default_rng(9)
        source=rng.uniform(0,20,(60,60)).astype('float32'); source[20:25,20:25]=np.nan;source[40:44,40:44]=0
        proxy=rng.uniform(.05,1,source.shape).astype('float32')
        ref=op.apply_fork_form_allocation(source,proxy,kernel=kernel,
                 denominator_epsilon_relative=1e-6,denominator_instability_threshold_relative=1e-4)
        gain=op.allocation_gain_from_components(ref.normalized_proxy,ref.convolved_proxy,
                 valid_mask=ref.geometric_support_fraction>=1-1e-6,denominator_epsilon=ref.denominator_epsilon)
        got,mean=allocate_nightly(source,gain,kernel,op._correlate_constant_zero,ref.geometric_support_fraction)
        np.testing.assert_allclose(got,ref.allocation,rtol=2e-6,atol=1e-6,equal_nan=True)
        np.testing.assert_equal(np.isfinite(got),ref.valid_output_mask)
        self.assertTrue(np.isnan(got[20:25,20:25]).all())
        empty,_=allocate_nightly(np.full(source.shape,np.nan),gain,kernel,op._correlate_constant_zero,ref.geometric_support_fraction)
        self.assertTrue(np.isnan(empty).all())

    @unittest.skipUnless((Path.home()/'research/ntl-psf-disaggregation/src').exists(),
                         'Original checkout required')
    def test_production_fft_kernel_and_cached_gain_match_upstream(self):
        m=upstream_modules(Path.home()/'research/ntl-psf-disaggregation');op=m['operator']
        k=m['kernels'].circular_mean_kernel(radius_m=500,resolution_m=10)
        rng=np.random.default_rng(44)
        y=rng.uniform(0,100,(256,256)).astype('float32')
        y[110:115,110:115]=np.nan
        h=rng.uniform(.05,1,y.shape).astype('float32')
        ref=op.apply_fork_form_allocation(y,h,kernel=k,
                denominator_epsilon_relative=1e-6,denominator_instability_threshold_relative=1e-4)
        normalized,norm=op.normalize_proxy_mean_one(h)
        corr=op._correlate_constant_zero
        geom=corr(np.ones_like(y),k.weights)
        denominator=corr(normalized,k.weights)/geom
        gain=op.allocation_gain_from_components(normalized,denominator,
                valid_mask=geom>=1-1e-6,denominator_epsilon=1e-6*norm.mean_after)
        allocated,uniform=allocate_nightly(y,gain,k,corr,geom)
        np.testing.assert_allclose(allocated,ref.allocation,rtol=2e-6,atol=1e-6,equal_nan=True)
        null=op.uniform_normalized_convolution_baseline(y,kernel=k,
                denominator_epsilon_relative=1e-6,denominator_instability_threshold_relative=1e-4)
        np.testing.assert_allclose(uniform,null.allocation,rtol=2e-6,atol=1e-6,equal_nan=True)

    @unittest.skipUnless((Path.home()/'research/ntl-psf-disaggregation/src').exists(),
                         'Original checkout required')
    def test_original_operator_is_not_cellwise_conserving(self):
        m=upstream_modules(Path.home()/'research/ntl-psf-disaggregation');op=m['operator']
        k=m['kernels'].circular_mean_kernel(radius_m=3,resolution_m=1)
        coarse=np.ones((9,9),dtype='float32');coarse[4,4]=10
        fine=np.repeat(np.repeat(coarse,10,axis=0),10,axis=1)
        proxy=np.ones_like(fine)
        # h=1 is the mandatory smoothing null and already disproves identity.
        result=op.apply_fork_form_allocation(fine,proxy,kernel=k,
                    denominator_epsilon_relative=1e-6,denominator_instability_threshold_relative=1e-4)
        center=result.allocation[40:50,40:50].mean()
        self.assertLess(center,10)
        self.assertGreater(center,1)

if __name__=='__main__': unittest.main()
