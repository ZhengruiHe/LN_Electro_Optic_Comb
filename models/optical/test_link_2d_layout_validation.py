"""实际芯层选择的安全测试；不依赖本机仿真结果，不启动求解器。"""
import unittest
from shapely.geometry import box
from shapely.ops import unary_union
from prepare_link_2d_layout_validation import identify_actual_pair


class ActualLayoutTests(unittest.TestCase):
    def setUp(self):
        self.pair = unary_union([box(-10, -.35, 10, .35), box(-10, 2.9, 10, 3.6)])
        self.roi = box(-9, -40, 9, 40)

    def test_exact_pair(self):
        result = identify_actual_pair(self.pair, self.pair, self.roi)
        self.assertTrue(result['actual_core_identified'])
        self.assertEqual(result['route_to_actual_max_shape_difference_um'], 0)

    def test_distant_third_route_is_explicitly_recorded(self):
        actual = unary_union([self.pair, box(-10, 22.9, 10, 23.6)])
        result = identify_actual_pair(actual, self.pair, self.roi)
        self.assertTrue(result['actual_core_identified'])
        self.assertAlmostEqual(result['excluded_other_route_min_edge_gap_um'], 19.3)

    def test_close_third_route_is_not_silently_ignored(self):
        actual = unary_union([self.pair, box(-10, 7, 10, 7.7)])
        self.assertFalse(identify_actual_pair(actual, self.pair, self.roi)['actual_core_identified'])

    def test_missing_waveguide_rejected(self):
        self.assertFalse(identify_actual_pair(box(-10, -.35, 10, .35), self.pair,
                                             self.roi)['actual_core_identified'])


if __name__ == '__main__':
    unittest.main()
