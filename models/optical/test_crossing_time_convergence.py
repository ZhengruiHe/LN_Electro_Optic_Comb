"""时间窗稳定和自动能量停止必须分别报告，不能由一个指标代替另一个。"""
import copy
import unittest

from check_crossing_time_convergence import compare_results, validate_baseline
from crossing_fdtd import MODEL_VERSION


def sample(status=1):
    return {"model_version": MODEL_VERSION, "spectral_validation_pass": True,
            "fdtd_status": status, "source": "north", "ports": {
                "north": {"S": {"real": .004, "imag": 0}, "power_db": -48.,
                          "actual_wavelength_nm": 1550.},
                "south": {"S": {"real": .967, "imag": 0}, "power_db": -.288,
                          "actual_wavelength_nm": 1550.}}}


class TimeTests(unittest.TestCase):
    def test_stable_but_no_energy_stop_is_not_pass(self):
        check = compare_results(sample(), sample())
        self.assertTrue(check["time_window_loss_stable"])
        self.assertFalse(check["time_check_pass"])

    def test_stable_and_energy_stop(self):
        self.assertTrue(compare_results(sample(), sample(2))["time_check_pass"])

    def test_energy_stop_without_stability_is_not_pass(self):
        new = copy.deepcopy(sample(2))
        new["ports"]["south"]["power_db"] -= .03
        self.assertFalse(compare_results(sample(), new)["time_check_pass"])

    def test_reject_old_frequency(self):
        old = sample()
        old["ports"]["south"]["actual_wavelength_nm"] = 1546.8387
        with self.assertRaises(ValueError):
            validate_baseline(old, {"time_fs": 2500., "source": "north", "wavelength_nm": 1550.}, 5000.)

    def test_reject_not_longer(self):
        with self.assertRaises(ValueError):
            validate_baseline(sample(), {"time_fs": 2500., "source": "north", "wavelength_nm": 1550.}, 2500.)


if __name__ == "__main__":
    unittest.main()
