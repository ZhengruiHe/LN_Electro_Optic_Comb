"""成对近邻弯的复模重构与局域TE0投影；不把超模编号直接当作波导编号。"""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def canonical(p,number):
    y=p['y'].ravel()*1e6
    e=p[f'E{number}'].squeeze()[:,1].astype(complex)
    h=p[f'H{number}'].squeeze()[:,2].astype(complex)
    power=.5*np.real(np.trapezoid(e*h.conj(),y))
    if power<0:h=-h;power=-power
    return y,e/np.sqrt(power),h/np.sqrt(power)


def overlap(y,e,h,em,hm):
    return .25*np.trapezoid(e*hm.conj()+em.conj()*h,y)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',type=Path,required=True);ap.add_argument('--reference',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    r=json.loads((a.case/'结果.json').read_text(encoding='utf-8'))
    geo=json.loads((a.case/'参数与几何.json').read_text(encoding='utf-8'));pair=geo['pair']
    state=json.loads((a.case/'状态.json').read_text(encoding='utf-8'))
    if state['status']!='component_completed':raise RuntimeError('联合场尚未完成')
    ds=r['ports']['near_pair']
    if not np.allclose(np.array(ds['lambda'])*1e9,1550,atol=1e-6,rtol=0):raise ValueError('实际波长不是1550nm')
    s=np.array([complex(v['real'],v['imag']) for v in ds['S'][0]])
    p=np.load(a.case/'near_pair_mode_profiles.npz')
    y,e,h=canonical(p,1);eo=np.zeros_like(e);ho=np.zeros_like(h)
    modal_branch=[]
    for i,amp in enumerate(s,1):
        parts=np.sum(abs(p[f'E{i}'].reshape(-1,3))**2,axis=0)
        planar=(parts[0]+parts[1])/sum(parts)
        modal_branch.append({'mode_number':i,'inplane_fraction':float(planar),'power':float(abs(amp)**2)})
        if planar<.9:continue
        _,ei,hi=canonical(p,i);eo+=amp*ei;ho+=amp*hi
    ref=np.load(a.reference/'input_mode_profiles.npz');yr,er,hr=canonical(ref,1)
    if not np.allclose(ref['lambda']*1e9,1550,atol=1e-6,rtol=0):
        raise ValueError('直参考模式频点不是1550nm')
    ref_state=json.loads((a.reference/'状态.json').read_text(encoding='utf-8'))
    ref_result=json.loads((a.reference/'结果.json').read_text(encoding='utf-8'))
    if ref_state['model_sha256']!=state['model_sha256'] or ref_result['fdtd_status']!=2:
        raise ValueError('直参考材料或能量停止状态不匹配')
    ref_s=ref_result['ports']['output']['S'][0][0]
    reference_power=abs(complex(ref_s['real'],ref_s['imag']))**2
    local=[]
    for center in pair['near_centers_um']:
        ei=np.interp(y-center,yr,er,left=0,right=0);hi=np.interp(y-center,yr,hr,left=0,right=0)
        power=.5*np.real(np.trapezoid(ei*hi.conj(),y))
        local.append((ei/np.sqrt(power),hi/np.sqrt(power)))
    gram=np.array([[overlap(y,ea,ha,eb,hb) for eb,hb in local] for ea,ha in local])
    gram=.5*(gram+gram.conj().T)
    vals,vecs=np.linalg.eigh(gram)
    if min(vals)<.9:raise RuntimeError('局域模式非正交程度过大，不能直接作两个近邻端口解释')
    projection=np.array([overlap(y,eo,ho,ei,hi) for ei,hi in local])
    # 弱耦合局域基底的对称正交化，保留原始Gram与投影供复核。
    a_local=(vecs@np.diag(vals**-.5)@vecs.conj().T)@projection
    power=abs(a_local)**2
    source=state['args']['source_port'];wanted=0 if source=='far_a' else 1
    def port_power(name):
        return sum(abs(complex(v['real'],v['imag']))**2 for v in r['ports'][name]['S'][0])
    reflection=port_power(source)
    backward_other=port_power('far_b' if source=='far_a' else 'far_a')
    field=np.load(a.case/'场分布.npz');xx=field['x'].ravel()*1e6;yy=field['y'].ravel()*1e6
    ey=field['E'].squeeze()[np.argmin(abs(xx)),:,1]
    sampled=np.interp(y,yy,ey,left=0,right=0)
    # 比较场形，不使用这里的绝对幅度标定。该监视器以3倍间隔采样。
    shape_similarity=float(abs(np.trapezoid(sampled.conj()*eo,y))**2/
        (np.trapezoid(abs(sampled)**2,y)*np.trapezoid(abs(eo)**2,y)))
    result={'scope':'等效2D近邻联合场；局域TE0弱耦合投影，不是全链路/三维验证',
            'pair':pair['pair'],'route_names':pair['route_names'],'source_port':source,
            'target_route':pair['route_names'][wanted],'wavelength_nm':1550.,'fdtd_status':r['fdtd_status'],
            'near_core_edge_gap_um':abs(pair['near_centers_um'][0]-pair['near_centers_um'][1])-.7,
            'power_by_local_TE0':power.tolist(),'tracked_near_port_mode_power':float(sum(abs(s)**2)),
            'local_Gram_abs':abs(gram).tolist(),'direct_field_vs_modal_reconstruction_similarity':shape_similarity,
            'local_target_insertion_loss_db':float(-10*np.log10(power[wanted])),
            'straight_reference_power':float(reference_power),
            'reference_normalized_target_insertion_loss_db':float(-10*np.log10(power[wanted]/reference_power)),
            'local_crosstalk_relative_db':float(10*np.log10(power[1-wanted]/power[wanted])),
            'tracked_reflection_power':float(reflection),
            'tracked_backward_other_route_power':float(backward_other),
            'unaccounted_power_in_selected_port_modes':float(1-sum(abs(s)**2)-reflection-backward_other),
            'power_balance_note':'未计入已选端口模的功率含辐射、未跟踪模式和数值误差，不直接称真实泄漏损耗',
            'source_gds_sha256':pair.get('source_gds_sha256'),
            'basis_checks_pass':bool(shape_similarity>.97 and abs(gram[0,1])<.05),
            'modal_branches':modal_branch,'full_chain_pass':False,
            'limitations':['单个入射路径，另一入射需另算','局域导模使用同网格直参考，并作弱重叠对称正交化',
                           '真实LN2映射及竖直偏振混合未由此模型验证']}
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    intensity=np.sum(abs(field['E'].squeeze())**2,axis=-1)
    fig,ax=plt.subplots(figsize=(8,7),layout='constrained')
    im=ax.pcolormesh(xx,yy,10*np.log10(np.maximum(intensity/intensity.max(),1e-7)).T,vmin=-45,vmax=0,shading='auto',cmap='turbo')
    for v in pair['core_polygons_um']:
        pts=np.array(v);ax.plot(pts[:,0],pts[:,1],color='white',lw=.4)
    labels={'left_lower':'左下接入','left_upper':'左上接入','right_upper':'右上接入','right_lower':'右下接入'}
    b=geo['bounds_um']
    ax.set(aspect='equal',xlim=(b[0],b[2]),ylim=(b[1],b[3]),xlabel='相对x（µm）',ylabel='相对y（µm）',
           title='实际成对接入弯：1550 nm等效2D联合场\n'+labels[pair['pair']]+'，'+pair['route_names'][wanted]+'激励')
    fig.colorbar(im,ax=ax,label='归一化 |E|²（dB）');fig.savefig(a.output_dir/'成对接入_二维联合传播场.png',dpi=160);plt.close(fig)
    with (a.output_dir/'近邻局域模式投影.json').open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
