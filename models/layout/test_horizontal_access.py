"""横向接入、独立周期与近邻几何回归；不运行电磁求解。"""
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_route_proximity import analyze
from route_horizontal_access import build
from route_geometry import check_route_network


class HorizontalRouteTests(unittest.TestCase):
    def test_cycle_polarity_guard(self):
        for p in [(3, 3, 4), (3.5, 3.5, 4), (3, 3.5, 3.5)]:
            with self.assertRaises(ValueError):
                build(periods=p)

    @unittest.skipUnless((ROOT / 'results/layout/_历史归档/迭代版本_20260912/YSJ合并15mm_功能候选_v10/合并版时延预算.json').exists(),
                         '需要本地预算输入')
    def test_horizontal_geometry_and_independent_cycles(self):
        r = build()
        self.assertEqual(r['loop_periods'], [3, 3.5, 4])
        self.assertEqual(r['electrode_length_um'], 15000)
        self.assertEqual(len(r['crossing_details']), 9)
        specs = []
        for c in r['crossing_details']:
            self.assertAlmostEqual(c['angle_deg'], 90, places=5)
            self.assertTrue(all(c['straight_245um_each_route'].values()))
            specs.append({'routes': c['routes'], 'center_um': [c['x_um'], c['y_um']], 'angle_deg': 90})
        self.assertTrue(check_route_network(r['routes'], specs)['topology_screen_pass'])
        for i, t in enumerate([300, 350, 400], 1):
            d = r['delays'][f'loop{i}']
            self.assertAlmostEqual(d['fixed_delay_ps'] + d['directional_delay_ps'], t, places=8)
        for c in r['shape_checks'].values():
            self.assertTrue(c['window_inside_block'])
            for key in ('YSJ_core_overlap_um2', 'SiN_removal_over_YSJ_SiN_um2', 'LN1_over_YSJ_LN2_etch_um2'):
                self.assertEqual(c[key], 0)
        report = analyze(r)
        self.assertTrue(report['ordinary_pairwise_clear'])
        self.assertTrue(report['no_close_long_parallel'])
        self.assertGreaterEqual(report['minimum_long_parallel_nominal_edge_gap_um'], 19.299)
        self.assertFalse(report['optical_crosstalk_verified'])
        # 固定端口附近仍有短程近邻，不能把其余路由通过写成全器件无耦合。
        self.assertLess(min(p['minimum_nominal_edge_gap_um'] for p in report['pairs']), 10)

    def test_parallel_proximity_is_detected(self):
        r = {'routes': {'a': [[0, 0], [1000, 0]], 'b': [[0, 3.25], [1000, 3.25]]},
             'crossing_details': []}
        report = analyze(r)
        self.assertFalse(report['ordinary_pairwise_clear'])
        self.assertFalse(report['no_close_long_parallel'])
        self.assertAlmostEqual(report['minimum_long_parallel_nominal_edge_gap_um'], 2.55)


if __name__ == '__main__':
    unittest.main()
