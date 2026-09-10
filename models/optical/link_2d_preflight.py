"""整链路降维前的MODE varFDTD直波导接口/各向异性校准检查。

只建立独立小工程，不修改任何版图或旧仿真；失败时保留工程和错误。
2.5D折叠到平面后不能自动当作完整矢量模式/偏振验证。
"""
from __future__ import annotations
import argparse
import json
import time
import traceback
from pathlib import Path

import numpy as np
import psutil
from shapely.geometry import box
from crossing_fdtd import add_layer, serial, assert_wavelength
from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi, zelmon_ln_indices


def save_json(path, value):
    path.write_text(json.dumps(serial(value), ensure_ascii=False, indent=2), encoding='utf-8')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--width-um', type=float, default=.7)
    p.add_argument('--axis', choices=['Y', 'Z'], default='Y')
    p.add_argument('--isotropic-control', action='store_true')
    p.add_argument('--mode', type=int, default=1)
    p.add_argument('--propagation-y', action='store_true', help='只旋转直波导，不旋转材料张量，用于核对同一平面模型的方向响应')
    p.add_argument('--run', action='store_true')
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=False)
    state = {'status': 'opening', 'controller_pid': psutil.Process().pid,
             'full_chain_pass': False, 'scope': '小直波导降维方法预检查，不是整链路结果',
             'args': {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}}
    def status(**kw):
        state.update(kw)
        state['updated_at'] = time.strftime('%Y-%m-%d %H:%M:%S')
        save_json(a.output_dir / '状态.json', state)
        print(state['status'], flush=True)
    status()
    cfg = load_config(DEFAULT_CONFIG)
    api = load_lumapi(cfg)
    try:
        with api.MODE(hide=True) as m:
            status(status='building', version=m.version())
            m.addvarfdtd()
            for key, value in {'x min': -8e-6, 'x max': 8e-6, 'y min': -8e-6, 'y max': 8e-6,
                               'z min': -1.5e-6, 'z max': 3.2e-6, 'simulation time': .8e-12,
                               'auto shutoff min': 1e-6, 'mesh accuracy': 2, 'polarization': 'E mode (TE)',
                               'x0': 0., 'y0': 0., 'can optimize mesh algorithm for extruded structures': False,
                               'bandwidth': 'narrowband'}.items():
                m.set(key, value)
            m.addrect()
            for key, value in {'name': 'SiO2_nominal', 'x span': 22e-6, 'y span': 22e-6,
                               'z min': -4e-6, 'z max': 2.1e-6,
                               'material': cfg['material_models']['sio2'],
                               'override mesh order from material database': True, 'mesh order': 3}.items():
                m.set(key, value)
            no, ne = [float(x[0]) for x in zelmon_ln_indices(np.array([1.55]))]
            # 光路始终沿仿真x，改变全局张量映射用于两个晶向校准。
            values = [no, ne, no] if a.axis == 'Y' else [ne, no, no]
            if a.isotropic_control:
                values = [ne, ne, ne]
            index = ';'.join(f'{x:.12g}' for x in values)
            slab = box(-3.6, -11, 3.6, 11) if a.propagation_y else box(-11, -3.6, 11, 3.6)
            core = box(-a.width_um/2, -11, a.width_um/2, 11) if a.propagation_y else box(-11, -a.width_um/2, 11, a.width_um/2)
            add_layer(m, slab, .7, .2, 8, 70, index, 'LN_residual_nominal7p2')
            add_layer(m, core, .9, .2, 8, 70, index, 'LN_rib')
            m.addmesh()
            for key, value in {'name': 'local_mesh', 'x span': 16e-6, 'y span': 16e-6,
                               'z min': .6e-6, 'z max': 1.2e-6, 'dx': 50e-9, 'dy': 50e-9, 'dz': 25e-9}.items():
                m.set(key, value)
            m.addmodesource()
            direction_axis, transverse = ('y', 'x') if a.propagation_y else ('x', 'y')
            for key, value in {'name': 'source', 'injection axis': direction_axis, 'direction': 'Forward',
                               direction_axis: -5e-6, transverse: 0., transverse+' span': 10e-6, 'mode selection': 'user select',
                               'override global source settings': False}.items():
                m.set(key, value)
            fc = 299792458 / 1.55e-6
            # MODE varFDTD源模按波长端点中点计算；此处必须以波长对称。
            # 宽源谱仅为短脉冲；所有保存的DFT结果均严格在1550nm。
            m.setglobalsource('set wavelength', True)
            m.setglobalsource('wavelength start', 1.5e-6)
            m.setglobalsource('wavelength stop', 1.6e-6)
            m.setglobalmonitor('use source limits', False)
            m.setglobalmonitor('frequency center', fc)
            m.setglobalmonitor('frequency span', 0.)
            m.setglobalmonitor('frequency points', 1)
            project = (a.output_dir / '直波导_varFDTD_降维预检查.lms').resolve()
            m.save(str(project))
            m.select('source')
            status(status='source_mode_calculation')
            updated = m.updatesourcemode(a.mode)
            source_mode = m.getresult('source', 'neff')
            assert_wavelength(source_mode, 1550., 'varFDTD source neff')
            source_profile = m.getresult('source', 'mode profile')
            assert_wavelength(source_profile, 1550., 'varFDTD source mode profile')
            np.savez_compressed(a.output_dir / 'source_mode_profile.npz',
                                **{k:v for k,v in source_profile.items() if isinstance(v,np.ndarray)})
            info = {'tensor_index_xyz': values, 'wavelength_nm': 1550., 'source_update_return': updated,
                    'source_mode_result': source_mode, 'source_wavelength_verified': True,
                    'available_results': {},
                    'assumptions': ['沿用400nm总LN/200nm脊/200nm有限7.2um残余平台、70度侧壁名义模型',
                                    '有限平台是旧仿真假设，不由当前GDS的20/1窗口直接推定',
                                    '窄带模型只用于1550nm模式，不用于重新估算ng']}
            for name in ['source', 'varFDTD']:
                try:
                    info['available_results'][name] = m.getresult(name)
                except Exception as exc:
                    info['available_results'][name] = {'unavailable': str(exc)}
            save_json(a.output_dir / '接口检查.json', info)
            m.addpower()
            for key, value in {'name': 'through', 'monitor type': '2D '+direction_axis.upper()+'-normal',
                               direction_axis: 5e-6, transverse: 0., transverse+' span': 10e-6}.items():
                m.set(key, value)
            m.addpower()
            for key, value in {'name': 'field_xy', 'monitor type': '2D Z-normal',
                               'x span': 14e-6, 'y span': 10e-6, 'z': 1e-6}.items():
                m.set(key, value)
            m.save(str(project))
            status(status='built', project=str(project))
            if a.run:
                previous = {}
                try:
                    for key, value in [('processes', '1'), ('threads', '4')]:
                        previous[key] = m.getresource('varFDTD', 1, key)
                        m.setresource('varFDTD', 1, key, value)
                    status(status='running')
                    start = time.time()
                    m.run()
                finally:
                    for key, value in previous.items():
                        m.setresource('varFDTD', 1, key, value)
                m.save(str(project))
                result = {'T': m.getresult('through', 'T'), 'solver_status': m.getresult('varFDTD', 'status'),
                          'elapsed_s': time.time()-start, 'full_chain_pass': False}
                assert_wavelength(result['T'], 1550., 'varFDTD through')
                save_json(a.output_dir / '结果.json', result)
                field = m.getresult('field_xy', 'E')
                assert_wavelength(field, 1550., 'varFDTD field')
                np.savez_compressed(a.output_dir / '场分布.npz', **{k:v for k,v in field.items() if isinstance(v,np.ndarray)})
                status(status='preflight_run_completed', elapsed_s=result['elapsed_s'])
    except Exception as exc:
        status(status='failed', error=str(exc), traceback=traceback.format_exc())
        raise


if __name__ == '__main__':
    main()
