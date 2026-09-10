"""在1550nm拟合面内各向异性的2D Hz等效模型；不冒充真实LN材料。

依据既有全矢量FDE的两晶向0.7um TE0、Y向1.33um TE0/TE1相位折射率。
解 d_y(eps_x^-1 d_y Hz)+(k0^2-beta^2/eps_y)Hz=0，固定全局晶轴。
上下厚度/70度侧壁通过FDE参考折合，不再给2D材料直接填写bulk LN折射率。
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
K0 = 2*np.pi/1.55


def solve_modes(params, width, axis='Y', dy=.02, count=5, n_clad=1.444):
    slab_x, slab_y, delta_x, delta_y = params
    core_x, core_y = slab_x+delta_x, slab_y+delta_y
    y = np.arange(-8, 8+dy/2, dy)
    def profile(core, slab):
        return np.where(abs(y) < width/2-1e-9, core**2,
                        np.where(abs(y) < 3.6-1e-9, slab**2, n_clad**2))
    eps_x, eps_y = profile(core_x, slab_x), profile(core_y, slab_y)
    if axis == 'Z':
        eps_x, eps_y = eps_y, eps_x
    # Hz未知量的端点为零，边界距核心足够远；这是直模实本征值，不计算辐射。
    inverse_edges = .5*(1/eps_x[:-1]+1/eps_x[1:])/dy**2
    ex = eps_y[1:-1]
    diagonal = (K0**2-inverse_edges[:-1]-inverse_edges[1:])*ex
    off = inverse_edges[1:-1]*np.sqrt(ex[:-1]*ex[1:])
    vals, vecs = eigh_tridiagonal(diagonal, off, select='i',
                                  select_range=(len(diagonal)-count, len(diagonal)-1))
    vals, vecs = vals[::-1], vecs[:, ::-1]
    neff = np.sqrt(np.maximum(vals, 0))/K0
    hz = np.zeros((len(y), count))
    hz[1:-1] = np.sqrt(ex)[:,None]*vecs
    ey = hz/eps_y[:,None]  # 横向电场除去每个模式不影响归一化的beta常数。
    ey /= np.sqrt(np.trapezoid(abs(ey)**2, y, axis=0))[None,:]
    return {'y_um': y, 'neff': neff, 'Hz': hz, 'E_transverse': ey}


def source_data():
    directional = ROOT/'results/optical/欧拉时延输入_0p7um方向群折射率_v1/0p7um_rib_两个晶向群折射率.json'
    active = ROOT/'results/optical/扫描结果/波导截面/pdk_1550_w1p33_e0p2_a70_ln2w7p2_fullsio2_strict.csv'
    dr = json.loads(directional.read_text(encoding='utf-8'))
    with active.open(encoding='utf-8') as f:
        rows = {r['mode_label']:r for r in csv.DictReader(f)}
    targets = np.array([dr['axes']['Y']['neff_1550'], dr['axes']['Z']['neff_1550'],
                        float(rows['TE0']['neff_real']), float(rows['TE1']['neff_real'])])
    sources = [{'path':str(p.relative_to(ROOT)), 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
               for p in (directional, active)]
    return targets, sources


def predict(params, dy=.02):
    a = solve_modes(params,.7,'Y',dy)['neff']
    b = solve_modes(params,.7,'Z',dy)['neff']
    c = solve_modes(params,1.33,'Y',dy)['neff']
    return np.array([a[0],b[0],c[0],c[1]])


def predict_routing(params):
    values=[]
    for axis in ('Y','Z'):
        r=solve_modes(params,.7,axis)
        y,e=r['y_um'],r['E_transverse'][:,0]
        values.extend([r['neff'][0],np.sqrt(np.trapezoid(y*y*abs(e)**2,y))])
    return np.array(values)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--routing-only',action='store_true',help='拟合0.7um两晶向的neff和场宽；有源/MUX复用矢量结果')
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    targets,sources=source_data()
    labels=['Y_0p7_TE0','Z_0p7_TE0','Y_1p33_TE0','Y_1p33_TE1']
    predictor=predict
    weights=np.ones(4)
    if a.routing_only:
        dr=json.loads((ROOT/sources[0]['path']).read_text(encoding='utf-8'))
        targets=np.array([v for axis in ('Y','Z') for v in (
            dr['axes'][axis]['neff_1550'],dr['axes'][axis]['samples'][1]['selected']['rms_x_um'])])
        labels=['Y_0p7_neff','Y_0p7_rms_um','Z_0p7_neff','Z_0p7_rms_um']
        predictor=predict_routing;weights=np.array([1.,.2,1.,.2])
    fit=least_squares(lambda p:(predictor(p)-targets)*weights,[1.69,1.64,.26,.25],
                      bounds=([1.445,1.445,.015,.015],[1.9,1.9,.5,.5]),
                      xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=300)
    sx,sy,dx,dy=fit.x
    predicted=predictor(fit.x)
    result={'scope':'拟合型1550nm各向异性2D模型；拟合吻合不是独立光学验证',
            'axes':{'x':'crystal_Y','y':'crystal_Z','z':'collapsed'},
            'polarization':'Hz，面内Ex/Ey；不预测真实竖直偏振杂化',
            'wavelength_nm':1550.,'parameters':fit.x.tolist(),
            'core_effective_index_xyz':[sx+dx,sy+dy,1.],
            'slab_effective_index_xyz':[sx,sy,1.],
            'cladding_index':1.444,'platform_top_width_um':7.2,
            'baseline_cross_section':'LN400nm/脊200nm/残余200nm，70度；有限7.2um平台为既有名义假设',
            'calibration_target':'0.7um两晶向的相位折射率和RMS场宽' if a.routing_only else '四个相位折射率',
            'usage':'仅无源路由/接口；有源段与MUX复用全矢量结果' if a.routing_only else '探索候选，尚未验证模式场形',
            'calibration':[{ 'metric':n,'target':float(t),'fitted':float(v),
                            'residual':float(v-t)} for n,t,v in zip(labels,targets,predicted)],
            'fit_success':bool(fit.success),'fit_max_neff_error':float(np.max(abs(predicted-targets))),
            'source_inputs':sources,'requires_independent_checks':['横向场形','网格','未参与拟合的波导宽度/方向'],
            'full_chain_pass':False}
    with (a.output_dir/'等效2D张量拟合.json').open('x',encoding='utf-8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    for width,axis in [(.7,'Y'),(.7,'Z'),(1.33,'Y'),(1.2,'Y'),(1.2,'Z')]:
        r=solve_modes(fit.x,width,axis)
        np.savez_compressed(a.output_dir/f'{axis}_w{width:g}_等效模式.npz',**r)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if not fit.success or result['fit_max_neff_error']>.005:
        raise RuntimeError('2D拟合未满足0.005相位折射率初始误差门槛，不得直接用于链路')


if __name__=='__main__':main()
