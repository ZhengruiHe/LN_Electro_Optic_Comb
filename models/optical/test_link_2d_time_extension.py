"""时间窗补算选择逻辑；不启动求解器。"""
import unittest
from run_link_2d_neighbors import select_time_extensions


class TimeExtensionTests(unittest.TestCase):
    def setUp(self):
        self.cases = [{'name': name, 'time_ps': 10} for name in ('a', 'b', 'c')]
        self.previous = {'status': 'neighbors_screening_completed', 'cases': self.cases,
                         'completed': [{'case': 'a', 'fdtd_status': 2}, {'case': 'b', 'fdtd_status': 1},
                                       {'case': 'c', 'fdtd_status': 2}]}

    def test_select_only_time_limited_case(self):
        selected = select_time_extensions(self.cases, self.previous, 16)
        self.assertEqual([x['name'] for x in selected], ['b'])
        self.assertEqual(selected[0]['time_ps'], 16)
        self.assertEqual(self.cases[1]['time_ps'], 10)

    def test_reject_live_batch(self):
        self.previous['status'] = 'running'
        with self.assertRaises(ValueError):
            select_time_extensions(self.cases, self.previous, 16)

    def test_reject_unchanged_window(self):
        with self.assertRaises(ValueError):
            select_time_extensions(self.cases, self.previous, 10)

    def test_no_duplicate_completed_cases(self):
        self.previous['completed'][1]['fdtd_status'] = 2
        with self.assertRaises(ValueError):
            select_time_extensions(self.cases, self.previous, 16)


if __name__ == '__main__':
    unittest.main()
