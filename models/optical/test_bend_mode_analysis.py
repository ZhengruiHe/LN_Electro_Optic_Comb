"""验证模式识别不依赖序号、处理两种端口轴、拒绝异常磁场和错误频点。"""
import unittest
import numpy as np
from bend_mode_analysis import mode_metrics


def synthetic(axis,broken_h=False):
    u=np.linspace(-5,5,101)*1e-6;z=np.linspace(-1,3,41)*1e-6
    field=np.exp(-(u[:,None]/.4e-6)**2-((z[None,:]-1e-6)/.3e-6)**2)
    profile={'lambda':np.array([1.55e-6]),'z':z}
    profile['y' if axis=='x' else 'x']=u
    trans=1 if axis=='x' else 0
    for i in [1,2]:
        E=np.zeros((101,41,3),complex);H=E.copy()
        # Mode1为TM；Mode2为TE0。
        E[:,:,2 if i==1 else trans]=field
        H[:,:,trans if i==1 else 2]=field*(1e-12 if broken_h else .005)
        profile['E'+str(i)]=E;profile['H'+str(i)]=H
    return profile,{'lambda':np.array([1.55e-6]),'neff':np.array([1.65,1.77])}


class ModeIdentityTests(unittest.TestCase):
    def test_axis_and_reordered_modes(self):
        for axis in ['x','y']:
            profile,neff=synthetic(axis)
            result=mode_metrics(profile,neff,axis)
            self.assertEqual(result['target_mode_number'],2)
            self.assertEqual(result['modes'][0]['label'],'TM_or_hybrid')

    def test_rejects_nonphysical_E_H_ratio(self):
        profile,neff=synthetic('x',broken_h=True)
        with self.assertRaises(ValueError):mode_metrics(profile,neff,'x')

    def test_rejects_wrong_wavelength(self):
        profile,neff=synthetic('x');profile['lambda']=np.array([1.551e-6])
        with self.assertRaises(ValueError):mode_metrics(profile,neff,'x')


if __name__=='__main__':unittest.main()
