"""重新规划的MUX S形拉锥几何测试。"""
import unittest
import numpy as np
from shapely.geometry import LineString, Polygon

from replan_mux_connections import minimum_radius, output_envelope_vertices, s_taper_vertices


class ReplanTests(unittest.TestCase):
    def test_ports_and_monotonic_center(self):
        for args in [(-200,0,0,-1.95,1.33,1.4),(750,950,-1.8,.15,1.1,.7),(750,950,1.45,3.4,.4,.7)]:
            x0,x1,c0,c1,w0,w1=args
            shape=Polygon(s_taper_vertices(*args))
            self.assertTrue(shape.is_valid)
            for x,c,w in ((x0,c0,w0),(x1,c1,w1)):
                cut=shape.intersection(LineString([(x,-10),(x,10)]))
                self.assertAlmostEqual(c,(cut.bounds[1]+cut.bounds[3])/2)
                self.assertAlmostEqual(w,cut.length)
            centers=[]
            for x in np.linspace(x0,x1,31):
                cut=shape.intersection(LineString([(x,-10),(x,10)]))
                centers.append((cut.bounds[1]+cut.bounds[3])/2)
            diff=np.diff(centers)
            self.assertTrue(np.all(diff<=1e-10) or np.all(diff>=-1e-10))

    def test_curvature_radius_above_pdk_minimum(self):
        for c0,c1 in ((0,-1.95),(-1.8,.15),(1.45,3.4)):
            self.assertGreater(minimum_radius(0,200,c0,c1),80)

    def test_output_envelope_parallel_then_matches_two_routes(self):
        for width in (7.2,17.2):
            shape=Polygon(output_envelope_vertices(width))
            self.assertTrue(shape.is_valid)
            for x,low,high in ((-200,-width/2,width/2),(750,-width/2,width/2),
                               (950,.15-width/2,3.4+width/2)):
                cut=shape.intersection(LineString([(x,-20),(x,20)]))
                self.assertAlmostEqual(cut.bounds[1],low)
                self.assertAlmostEqual(cut.bounds[3],high)


if __name__=='__main__':
    unittest.main()
