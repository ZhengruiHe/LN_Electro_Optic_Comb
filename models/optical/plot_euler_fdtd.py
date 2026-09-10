"""从真实FDTD结果绘制欧拉弯传播场和输出模式功率；不启动求解器。"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from bend_mode_analysis import analyze_folder


def plot(folder,output):
    folder=Path(folder);output=Path(output)
    if output.exists():raise FileExistsError('图文件已存在，拒绝覆盖')
    analysis=analyze_folder(folder)
    info=json.loads((folder/'参数与几何.json').read_text(encoding='utf-8'))
    state=json.loads((folder/'状态.json').read_text(encoding='utf-8'))
    with np.load(folder/'场分布.npz') as data:
        x=data['x'].ravel()*1e6;y=data['y'].ravel()*1e6
        E=data['E'].squeeze()
        if E.shape!=(len(x),len(y),3):raise ValueError('场数据维度错误')
        intensity=np.sum(abs(E)**2,axis=-1)
        if intensity.max()<=0:raise ValueError('场强为零')
        intensity/=intensity.max()
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                         'axes.unicode_minus':False})
    fig,(ax,bar)=plt.subplots(1,2,figsize=(12,5),gridspec_kw={'width_ratios':[1.25,1.]})
    image=ax.pcolormesh(x,y,intensity.T,shading='auto',cmap='inferno',vmin=0,vmax=1,rasterized=True)
    for key,color in [('platform_polygon_um','#00bcd4'),('core_polygon_um','#ffffff')]:
        pts=np.asarray(info[key]);ax.plot(pts[:,0],pts[:,1],color=color,linewidth=.65,alpha=.8)
    ax.set(aspect='equal',xlabel='全局x / 晶体Y (µm)',ylabel='全局y / 晶体Z (µm)',
           title='1550 nm真实几何传播场，z=1.0 µm')
    fig.colorbar(image,ax=ax,label='切片归一化 |E|²（不代表透射率）',fraction=.045,pad=.03)
    powers=np.asarray(state['output_power_by_mode']).ravel()
    db=10*np.log10(np.maximum(powers,1e-12))
    colors=['#0072b2' if i+1==analysis['output_TE0_number'] else '#d55e00' for i in range(len(powers))]
    bar.bar(np.arange(1,len(powers)+1),np.maximum(db,-80)+80,bottom=-80,align='center',color=colors)
    bar.scatter(np.arange(1,len(powers)+1),np.maximum(db,-80),c=colors,s=20,zorder=3)
    bar.set(ylim=(-80,0.1),xticks=np.arange(1,len(powers)+1),xlabel='输出端口模式编号',
            ylabel='模式功率 / 输入功率 (dB)',title='蓝色：场分布确认的TE0；橙色：其他模式')
    bar.grid(axis='y',alpha=.2)
    fig.suptitle(f"TE0原始插损 {analysis['raw_TE0_insertion_loss_db']:.4f} dB；"
                 f"停止状态 {analysis['fdtd_status']}；单次网格结果待收敛比较")
    fig.tight_layout();fig.savefig(output,dpi=180);plt.close(fig)
    return str(output.resolve())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(plot(a.folder,a.output))
