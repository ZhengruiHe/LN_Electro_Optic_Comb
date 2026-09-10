"""MUX包络平滑偏移的端点、宽度及对称性测试。"""
import unittest
import numpy as np
from shapely.geometry import LineString, Polygon

from align_mux_envelopes_layout import shifted_envelope_vertices


class EnvelopeTests(unittest.TestCase):
    def test_end_centers_and_widths(self):
        for width in (7.2,17.2):
            shape=Polygon(shifted_envelope_vertices(width))
            self.assertTrue(shape.is_valid)
            for x,center in ((-200,-1.95),(0,0),(950,0)):
                cut=shape.intersection(LineString([(x,-20),(x,20)]))
                self.assertAlmostEqual(cut.length,width)
                self.assertAlmostEqual((cut.bounds[1]+cut.bounds[3])/2,center)

    def test_transition_center_is_monotonic(self):
        shape=Polygon(shifted_envelope_vertices(7.2))
        centers=[]
        for x in np.linspace(-200,0,41):
            cut=shape.intersection(LineString([(x,-20),(x,20)]))
            centers.append((cut.bounds[1]+cut.bounds[3])/2)
        self.assertTrue(np.all(np.diff(centers)>=-1e-12))
        self.assertAlmostEqual(centers[0],-1.95)
        self.assertAlmostEqual(centers[-1],0)


if __name__=='__main__':
    unittest.main()
