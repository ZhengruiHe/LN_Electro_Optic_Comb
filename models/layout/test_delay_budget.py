"""欧拉时延方向积分与长度反解测试。"""
import unittest

from euler_routes import euler_return_route

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from delay_budget import route_delay, solve_length


class DelayTests(unittest.TestCase):
    def test_axis_aligned_delay(self):
        row=route_delay([[0,0],[100,0],[100,50]],2.2,2.3)
        self.assertAlmostEqual(row['crystal_Y_weighted_length_um'],100)
        self.assertAlmostEqual(row['crystal_Z_weighted_length_um'],50)

    def test_length_solver_closes_target(self):
        builder=lambda length:euler_return_route(x_right=18150,x_left=850,y_start=29.85,y_end=33.4,
            target_length=length,outward_sign=1,start_radius=80)
        length,points,bends,row,history=solve_length(builder,23447.4,172.0,2.2,2.3)
        self.assertAlmostEqual(row['directional_delay_ps'],172.0,places=8)
        self.assertGreaterEqual(bends['minimum_radius_um'],80)
        self.assertTrue(history)


if __name__=='__main__':
    unittest.main()
