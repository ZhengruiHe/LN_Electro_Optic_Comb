"""汇总已完成的2D阶段结果与独立截面比较；不把基础组件称整链路通过。"""
import argparse,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from link_2d_effective_model import solve_modes


def numbers(value):
    return np.array([complex(v['real'],v['imag']) for v in value[0]])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    model=json.loads((a.root/'无源各向异性拟合_v2/等效2D张量拟合.json').read_text(encoding='utf-8'))
    reference=json.loads((a.root/'独立FDE参考_w1p2/状态.json').read_text(encoding='utf-8'))
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    calibration=[]
    for col,axis in enumerate(('Y','Z')):
        modes=[v for v in reference['results'][axis]['modes'] if v['central_energy_fraction']>.7]
        reduced=solve_modes(model['parameters'],1.2,axis)
        for order,row in enumerate(modes[:2]):
            original=np.load(a.root/'独立FDE参考_w1p2'/f"{axis}_mode{row['solver_mode']}_全矢量场.npz")
            x=original['x'].ravel()*1e6;z=original['y'].ravel()*1e6
            field=sum(abs(original[k].squeeze())**2 for k in ['Ex','Ey','Ez'])
            intensity=np.trapezoid(field,z,axis=1)
            intensity/=np.trapezoid(intensity,x)
            yy=reduced['y_um'];ii=abs(reduced['E_transverse'][:,order])**2
            interp=np.interp(yy,x,intensity,left=0,right=0);interp/=np.trapezoid(interp,yy)
            similarity=float(np.trapezoid(np.sqrt(ii*interp),yy)**2)
            row2={'axis':axis,'width_um':1.2,'mode':'TE'+str(order),'FDE_neff':row['neff'],
                  'reduced_neff':float(reduced['neff'][order]),'intensity_profile_similarity':similarity,
                  'metric_note':'归一化横向强度轮廓相似度，不是完整矢量复场重叠或实测串扰'}
            row2['neff_relative_error_percent']=100*(row2['reduced_neff']/row2['FDE_neff']-1)
            calibration.append(row2)
            ax=axes[order,col];ax.plot(x,intensity,label='真实厚度/侧壁：矢量FDE');ax.plot(yy,ii,'--',label='等效2D模型')
            ax.set(xlim=(-2.5,2.5),xlabel='横向（µm）',ylabel='归一化横向强度',
                   title=f"晶体{axis}向 / TE{order}：轮廓相似度 {100*similarity:.2f}%")
            ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('1.2 µm独立截面检查：该宽度未参与0.7 µm模型拟合',fontsize=13)
    fig.savefig(a.output_dir/'等效模型_独立截面比较.png',dpi=160);plt.close(fig)
    completed=[]
    for directory in ('二维组件','基础组件批次_01','近邻联合仿真'):
        for p in sorted((a.root/directory).glob('*/结果.json')):
            r=json.loads(p.read_text(encoding='utf-8'))
            if 'output' not in r['ports']:
                completed.append({'case':p.parent.name,'result':str(p),'fdtd_status':r['fdtd_status'],
                                  'modal_projection_pending':True});continue
            for ds in r['ports'].values():
                if not np.allclose(np.array(ds['lambda'])*1e9,1550,atol=1e-6,rtol=0):raise ValueError('非1550nm数据')
            out=abs(numbers(r['ports']['output']['S']))**2
            back=abs(numbers(r['ports']['input']['S']))**2
            unwanted=float(sum(out[1:]))
            completed.append({'case':p.parent.name,'result':str(p),'fdtd_status':r['fdtd_status'],
                              'raw_target_power':float(out[0]),'raw_insertion_loss_db':float(-10*np.log10(out[0])),
                              'tracked_unwanted_output_power':unwanted,'tracked_output_mode_count':len(out),
                              'tracked_mode_suppression_db':float(10*np.log10(out[0]/max(unwanted,1e-300))),
                              'total_tracked_reflection':float(sum(back)),
                              'raw_missing_power':float(1-sum(out)-sum(back)),
                              'full_chain_pass':False})
    bend=a.root/'二维组件/R80_90度Y_100nm'
    if (bend/'场分布.npz').exists():
        d=np.load(bend/'场分布.npz');x=d['x'].ravel()*1e6;y=d['y'].ravel()*1e6
        intensity=np.sum(abs(d['E'].squeeze())**2,axis=-1)
        geo=json.loads((bend/'参数与几何.json').read_text(encoding='utf-8'))
        fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
        im=ax.pcolormesh(x,y,10*np.log10(np.maximum(intensity/intensity.max(),1e-6)).T,
                         vmin=-40,vmax=0,cmap='turbo',shading='auto')
        polys=geo.get('core_polygons_um',[geo.get('core_polygon_um')])
        for poly in polys:
            v=np.array(poly);ax.plot(v[:,0],v[:,1],color='white',lw=.4)
        ax.set(aspect='equal',xlabel='晶体Y方向（µm）',ylabel='晶体Z方向（µm）',
               title='R80欧拉弯：等效2D传播场，1550 nm\n不是三维弯曲或完整链路验证')
        fig.colorbar(im,ax=ax,label='归一化 |E|²（dB）')
        fig.savefig(a.output_dir/'R80欧拉弯_二维传播场.png',dpi=160);plt.close(fig)
    result={'status':'阶段结果；基础组件和近邻联合仿真仍在推进','wavelength_nm':1550.,
            'independent_cross_section':calibration,'completed_components':completed,
            'blackboxes':'9只交叉器按用户确认理想直通占位，不计其真实损耗/串扰',
            'limitations':['仅面内等效模式，不包含真实竖直TE/TM杂化',
                           '当前数据主要为各组件目标输入的一列复S，不是整链路完整多模散射矩阵',
                           '有限7.2um残余平台沿用旧名义仿真假设，未由未绘制的21层得到确认',
                           '不能用孤立弯通过替代四组近邻联合场和完整级联',
                           '原始插损需要直参考、网格及级联后再解释，负的小数值可能是数值归一化误差'],
            'full_chain_pass':False}
    with (a.output_dir/'阶段汇总.json').open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
