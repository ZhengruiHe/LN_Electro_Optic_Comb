"""二维链路方法的方程与几何单元测试，不启动求解器。"""
import unittest
import numpy as np
from link_2d_effective_model import solve_modes
from link_2d_component import geometry


class ReducedModelTests(unittest.TestCase):
    def test_isotropic_axes_agree(self):
        p=[1.65,1.65,.25,.25]
        np.testing.assert_allclose(solve_modes(p,.7,'Y')['neff'],solve_modes(p,.7,'Z')['neff'],atol=1e-10)

    def test_anisotropic_axes_remain_distinct(self):
        p=[1.7082380311642786,1.6644971880146087,.22380427792781757,.2010967256896744]
        n1=solve_modes(p,.7,'Y')['neff'][0];n2=solve_modes(p,.7,'Z')['neff'][0]
        self.assertGreater(n2-n1,.05)

    def test_leads_extend_beyond_ports_and_boundaries(self):
        for kind in ('reference','taper','bend90','bend180','s_bend'):
            for axis in ('Y','Z'):
                core,slab,pts,ports,bounds=geometry(kind,80,axis)
                self.assertTrue(core.is_valid)
                for normal,position,_,direction in ports.values():
                    i=0 if normal=='x' else 1
                    boundary=bounds[i] if direction=='Forward' else bounds[i+2]
                    self.assertAlmostEqual(abs(position-boundary),4.)
                    # 中心线在该侧至少越过边界2um，避免端面在PML前反射。
                    extent=pts[:,i].min() if direction=='Forward' else pts[:,i].max()
                    self.assertGreaterEqual(abs(extent-position),6.-1e-6)


if __name__=='__main__':unittest.main()
