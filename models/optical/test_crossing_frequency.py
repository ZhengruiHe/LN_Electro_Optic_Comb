"""避免把波长标签、端口实际频点和场监视器实际频点混为一谈。"""
import unittest
import numpy as np
from crossing_fdtd import assert_wavelength, configure_exact_frequency, C0


class FakeFrequencyAPI:
    def __init__(self):
        self.source, self.monitor, self.named = {}, {}, {}
    def setglobalsource(self,k,v): self.source[k]=v
    def setglobalmonitor(self,k,v): self.monitor[k]=v
    def setnamed(self,o,k,v): self.named[(o,k)]=v
    def getglobalsource(self,k):
        return (self.source['frequency start']+self.source['frequency stop'])/2
    def getglobalmonitor(self,k): return self.monitor[k]


class FrequencyTests(unittest.TestCase):
    def test_exact_frequency(self):
        f=FakeFrequencyAPI()
        configure_exact_frequency(f,1550.)
        self.assertAlmostEqual(C0/f.getglobalsource('center frequency')*1e9,1550.)
        self.assertEqual(f.monitor['frequency span'],0)
        self.assertEqual(f.named[('FDTD::ports','monitor frequency points')],1)
    def test_reject_old_port_wavelength(self):
        with self.assertRaisesRegex(ValueError,'实际波长'):
            assert_wavelength({'lambda':np.array([1.5468387096774196e-6])},1550,'port')
    def test_reject_old_field_wavelength(self):
        with self.assertRaises(ValueError):
            assert_wavelength({'lambda':[1.5354838709677424e-6]},1550,'field')
    def test_exact_sample(self):
        self.assertAlmostEqual(assert_wavelength({'lambda':[[1.55e-6]]},1550,'S'),1550)
    def test_multiple_samples_not_silently_first(self):
        with self.assertRaises(ValueError):
            assert_wavelength({'lambda':[1.55e-6,1.551e-6]},1550,'S')


if __name__=='__main__': unittest.main()
