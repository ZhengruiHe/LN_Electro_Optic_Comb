"""部分欧拉弯和三条四程回路的几何闭合测试。"""
import math
import unittest

from euler_routes import build_euler_routes, partial_euler_local
from route_geometry import check_route_network


TARGETS=[23447.44797283442,30258.420481281406,23447.44797283442]


class EulerTests(unittest.TestCase):
    def test_bend_angle_and_minimum_radius(self):
        for angle in (math.pi/2,-math.pi/2,math.pi,-math.pi):
            points,info=partial_euler_local(angle,80)
            self.assertAlmostEqual(info['end_heading_rad'],angle,places=8)
            self.assertAlmostEqual(info['maximum_curvature_per_um'],1/80)
            self.assertGreater(len(points),100)

    def test_euler_routes_close_and_keep_length(self):
        routes,info=build_euler_routes(TARGETS)
        expected={"loop1":([18150,29.85],[850,33.4]),"loop2":([18150,26.6],[850,-30.15]),
                  "loop3":([18150,-29.85],[850,-33.4]),"input":([0,250],[850,30.15])}
        for name,(start,end) in expected.items():
            self.assertAlmostEqual(routes[name][0][0],start[0],places=5)
            self.assertAlmostEqual(routes[name][0][1],start[1],places=5)
            self.assertAlmostEqual(routes[name][-1][0],end[0],places=5)
            self.assertAlmostEqual(routes[name][-1][1],end[1],places=5)
        for name,target in zip(("loop1","loop2","loop3"),TARGETS):
            self.assertAlmostEqual(info[name]['actual_centerline_length_um'],target,places=2)
            self.assertGreaterEqual(info[name]['minimum_radius_um'],80)
        self.assertEqual(info['total_euler_bends'],19)

    def test_three_delay_loops_do_not_cross(self):
        routes,_=build_euler_routes(TARGETS)
        report=check_route_network({name:routes[name] for name in ('loop1','loop2','loop3')})
        self.assertTrue(report['topology_screen_pass'],report)

    def test_no_long_diagonal_straight_connector(self):
        routes,_=build_euler_routes(TARGETS)
        for name,points in routes.items():
            for first,second in zip(points,points[1:]):
                dx=second[0]-first[0]
                dy=second[1]-first[1]
                length=math.hypot(dx,dy)
                if length>2.0:
                    self.assertTrue(abs(dx)<1e-9 or abs(dy)<1e-9,
                                    f'{name}出现长斜直线：{first}->{second}')

    def test_input_route_and_reserved_crossing_are_valid(self):
        routes,_=build_euler_routes(TARGETS)
        report=check_route_network(routes,[{"routes":["input","loop2"],"center_um":[300,250],"angle_deg":90}])
        self.assertTrue(report['topology_screen_pass'],report)


if __name__=='__main__':
    unittest.main()
