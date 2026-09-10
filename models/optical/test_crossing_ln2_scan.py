"""LN2扫描不改端口或LN1参数，也不把初筛排名当成收敛通过。"""
import unittest

from scan_crossing_ln2 import make_plan, summarize


class LN2ScanTests(unittest.TestCase):
    def test_plan_contains_both_axes_and_fixed_references(self):
        plan = make_plan([3.6,10.8,7.2],2500.)
        self.assertEqual(len(plan),8)
        for row in plan:
            a = row["params"]
            self.assertEqual((a["width_x_um"],a["width_y_um"],a["half_length_um"]),(2.5,2.5,16.))
            self.assertEqual(a["shutoff"],1e-7)
            if row["variant"] == "straight":
                self.assertEqual(a["slab_center_width_um"],7.2)

    def test_requires_original_baseline(self):
        with self.assertRaises(ValueError):
            make_plan([3.6,10.8],2500.)

    def test_invalid_time_rejected(self):
        with self.assertRaises(ValueError):
            make_plan([7.2],float("nan"))

    def test_status_one_never_becomes_a_pass(self):
        rows = []
        for source in ("north","west"):
            for variant, trans in (("crossing",.96),("straight",1.)):
                rows.append({"source":source,"variant":variant,"slab_center_width_um":7.2,
                             "transmission":trans,"insertion_loss_db":.18 if trans<1 else 0,
                             "reflection_db":-40.,"max_crosstalk_db":-50.,"energy_shutoff_reached":False})
        summary = summarize(rows)
        self.assertEqual(len(summary["preliminary_ranking"]),1)
        self.assertFalse(summary["preliminary_ranking"][0]["all_energy_shutoffs_reached"])
        self.assertFalse(summary["confirmed_nominal_pass"])


if __name__ == "__main__":
    unittest.main()
