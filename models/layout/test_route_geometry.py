"""运行：python -m unittest discover -s models/layout -p test_route_geometry.py。"""

import unittest

from route_geometry import check_route_network


class RouteGeometryTests(unittest.TestCase):
    def test_input_must_not_be_omitted(self):
        routes = {"loop1": [[0, 0], [2, 0]], "input": [[1, -1], [1, 1]]}
        report = check_route_network(routes)
        self.assertFalse(report["topology_screen_pass"])
        self.assertEqual(len(report["unexpected_intersections"]), 1)

    def test_explicit_orthogonal_crossing(self):
        routes = {"input": [[-1, 0], [1, 0]], "loop2": [[0, -1], [0, 1]]}
        reserved = [{"routes": ["input", "loop2"], "center_um": [0, 0]}]
        report = check_route_network(routes, reserved)
        self.assertTrue(report["topology_screen_pass"])
        self.assertFalse(report["optical_function_pass"])

    def test_wrong_crossing_angle_rejected(self):
        routes = {"input": [[-1, 0], [1, 0]], "loop2": [[-1, -1], [1, 1]]}
        reserved = [{"routes": ["input", "loop2"], "center_um": [0, 0]}]
        self.assertFalse(check_route_network(routes, reserved)["topology_screen_pass"])

    def test_reservation_does_not_whitelist_another_crossing(self):
        routes = {"input": [[-1, 0], [2, 0]], "loop2": [[0, -1], [0, 1]],
                  "output": [[1, -1], [1, 1]]}
        reserved = [{"routes": ["input", "loop2"], "center_um": [0, 0]}]
        self.assertFalse(check_route_network(routes, reserved)["topology_screen_pass"])

    def test_overlap_and_self_crossing_rejected(self):
        routes = {"input": [[0, 0], [2, 0]], "loop2": [[0.5, 0], [1.5, 0]]}
        self.assertFalse(check_route_network(routes)["topology_screen_pass"])
        routes = {"loop1": [[0, 0], [1, 1], [0, 1], [1, 0]]}
        self.assertFalse(check_route_network(routes)["topology_screen_pass"])

    def test_missing_reservation_rejected(self):
        routes = {"input": [[0, 0], [1, 0]], "loop2": [[0, 2], [1, 2]]}
        reserved = [{"routes": ["input", "loop2"], "center_um": [0, 0]}]
        self.assertFalse(check_route_network(routes, reserved)["topology_screen_pass"])


if __name__ == "__main__":
    unittest.main()
