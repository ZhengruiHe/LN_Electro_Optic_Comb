"""PDK前置门槛的无许可证单元测试。"""
import copy
import unittest
from types import SimpleNamespace

from shapely.geometry import LineString, box
from shapely.ops import unary_union

from crossing_fdtd import geometry
from crossing_pdk import validate_nominal
from mode_pdk_sweep import load_config, DEFAULT_CONFIG


class CrossingPDKTests(unittest.TestCase):
    def setUp(self):
        self.args = SimpleNamespace(width_x_um=2., width_y_um=2., port_width_um=.7,
                                    half_length_um=12., reference_axis=None)
        self.cfg = load_config(DEFAULT_CONFIG)

    def check(self):
        core, slab = geometry(self.args, self.args.half_length_um+6)
        return validate_nominal(self.args, self.cfg, core, slab)

    def test_all_planned_geometries(self):
        for w in (1.5, 2., 2.5, 3.):
            for length in (8., 12., 16.):
                self.args.width_x_um = self.args.width_y_um = w
                self.args.half_length_um = length
                self.assertTrue(self.check()["parameter_precheck_pass"])

    def test_sidewall_bottom_width(self):
        self.assertAlmostEqual(self.check()["port_rib_bottom_width_um"], .8455880937, places=8)

    def test_port_below_manual_limit(self):
        self.args.port_width_um = .29
        with self.assertRaisesRegex(ValueError, "最小线宽"):
            self.check()

    def test_core_bottom_outside_platform(self):
        self.args.width_x_um = self.args.width_y_um = 10.
        with self.assertRaisesRegex(ValueError, "超出LN2"):
            self.check()

    def test_config_change_is_not_ignored(self):
        for section, key, value in (("stack_um", "ln", .3),
                                    ("waveguide_um", "sidewall_angle_deg_from_horizontal", 60),
                                    ("waveguide_um", "active_region_sin_fully_removed", False)):
            original = copy.deepcopy(self.cfg)
            self.cfg[section][key] = value
            with self.assertRaisesRegex(ValueError, "不一致"):
                self.check()
            self.cfg = original

    def test_straight_reference(self):
        for axis in ("x", "y"):
            self.args.reference_axis = axis
            self.assertTrue(self.check()["parameter_precheck_pass"])

    def test_default_platform_matches_historical_geometry(self):
        _, slab = geometry(self.args, 18)
        expected = unary_union([box(-18,-3.6,18,3.6), box(-3.6,-18,3.6,18)])
        self.assertTrue(slab.equals(expected))

    def test_local_platform_width_preserves_core_and_port(self):
        self.args.width_x_um = self.args.width_y_um = 2.5
        self.args.half_length_um = 16.
        core_base, _ = geometry(self.args, 22)
        for width in (3.6,7.2,10.8):
            self.args.slab_center_width_um = width
            core, slab = geometry(self.args, 22)
            self.assertTrue(core.equals(core_base))
            for line in (LineString([(18,-10),(18,10)]), LineString([(-10,18),(10,18)])):
                self.assertAlmostEqual(slab.intersection(line).length, 7.2)
            self.assertTrue(self.check()["parameter_precheck_pass"])

    def test_reject_unsupported_rib_bottom(self):
        self.args.width_x_um = self.args.width_y_um = 2.5
        self.args.slab_center_width_um = 1.
        with self.assertRaisesRegex(ValueError, "超出LN2"):
            self.check()

    def test_reject_nan_platform(self):
        self.args.slab_center_width_um = float("nan")
        core, slab = geometry(SimpleNamespace(**{**vars(self.args), "slab_center_width_um":7.2}),18)
        with self.assertRaisesRegex(ValueError, "有限正数"):
            validate_nominal(self.args,self.cfg,core,slab)

    def test_reference_must_keep_original_platform(self):
        self.args.reference_axis = "x"
        self.args.slab_center_width_um = 3.6
        with self.assertRaisesRegex(ValueError, "直波导基准"):
            self.check()


if __name__ == "__main__":
    unittest.main()
