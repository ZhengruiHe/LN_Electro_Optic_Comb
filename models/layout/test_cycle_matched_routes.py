"""周期布线路由的回归检查；缺少本地求解结果时跳过数据依赖测试。"""
import unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from route_cycle_matched import build
from route_geometry import check_route_network

class CycleRouteTests(unittest.TestCase):
    def test_polarity_period_constraint(self):
        with self.assertRaises(ValueError):build(periods=(3.,3.,3.))

    def test_aux_tail_length_guard(self):
        with self.assertRaises(ValueError):build(active=20.)

    @unittest.skipUnless((ROOT/'results/layout/_历史归档/迭代版本_20260912/YSJ合并15mm_功能候选_v10/合并版时延预算.json').exists(),'需要本地预算输入')
    def test_complete_nominal_routes(self):
        r=build()
        self.assertEqual(len(r['routes']),5)
        self.assertEqual(len(r['crossing_details']),7)
        self.assertEqual(r['topology']['self_crossing_routes'],[])
        specs=[]
        for c in r['crossing_details']:
            self.assertAlmostEqual(c['angle_deg'],90.,places=5)
            self.assertTrue(all(c['straight_245um_each_route'].values()))
            specs.append({'routes':c['routes'],'center_um':[c['x_um'],c['y_um']],'angle_deg':90.})
        self.assertTrue(check_route_network(r['routes'],specs)['topology_screen_pass'])
        for i,t in enumerate([300.,350.,300.],1):
            d=r['delays'][f'loop{i}']
            self.assertAlmostEqual(d['fixed_delay_ps']+d['directional_delay_ps'],t,places=8)
        for c in r['shape_checks'].values():
            self.assertTrue(c['window_inside_block'])
            self.assertEqual(c['YSJ_core_overlap_um2'],0.)
            self.assertEqual(c['SiN_removal_over_YSJ_SiN_um2'],0.)
            self.assertEqual(c['LN1_over_YSJ_LN2_etch_um2'],0.)

if __name__=='__main__':unittest.main()
