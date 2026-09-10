"""黑盒外接锥形的端宽、平滑轨迹与空间限制测试。"""
import unittest
import numpy as np
from shapely.geometry import Polygon, LineString
from integrate_blackbox_crossing_layout import taper_vertices


class TaperTests(unittest.TestCase):
    def test_end_widths(self):
        shape = Polygon(taper_vertices(100.))
        self.assertTrue(shape.is_valid)
        for x,width in ((0,.7),(100,1.2)):
            self.assertAlmostEqual(shape.intersection(LineString([(x,-2),(x,2)])).length,width)

    def test_monotonic_and_symmetric(self):
        pts = np.asarray(taper_vertices(100.))
        widths = -2*pts[:201,1]
        self.assertTrue(np.all(np.diff(widths)>=0))
        np.testing.assert_allclose(pts[:201,0],pts[201:][::-1,0])
        np.testing.assert_allclose(pts[:201,1],-pts[201:][::-1,1])

    def test_reject_length_invading_bend(self):
        for length in (-1,0,125,float('nan')):
            with self.assertRaises(ValueError):
                taper_vertices(length)


if __name__ == '__main__':
    unittest.main()
