"""从保存的实际GDS回读几何生成MUX位置和窗口接头诊断图，不修改GDS。"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PlotPolygon
from shapely.geometry import box, LineString
from search_fixed_block_placement import geometry


def paint(ax, shape, color, edge, origin=(0.,0.), alpha=1.):
    for g in shape.geoms if hasattr(shape, 'geoms') else [shape]:
        if g.is_empty or g.geom_type!='Polygon':continue
        coords=[(x-origin[0],y-origin[1]) for x,y in g.exterior.coords]
        ax.add_patch(PlotPolygon(coords,facecolor=color,edgecolor=edge,lw=.7,alpha=alpha))
        for hole in g.interiors:
            ax.add_patch(PlotPolygon([(x-origin[0],y-origin[1]) for x,y in hole.coords],
                                    facecolor='white',edgecolor=edge,lw=.6))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--directory',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    data=json.loads((a.directory/'合并回读_四程几何.json').read_text(encoding='utf-8'))
    geoms={key:geometry(data,key) for key in [(20,0),(20,1),(10,2),(42,0)]}
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    fig=plt.figure(figsize=(13,8.5),layout='constrained')
    grid=fig.add_gridspec(2,2,height_ratios=[1.2,1])
    ax=fig.add_subplot(grid[0,:]);region=box(-10780,1570,-9490,1675)
    paint(ax,geoms[(42,0)].intersection(region),'#edb6b0','#c7645d',alpha=.35)
    ax.axvspan(-10470,-9720,color='#a8ceef',alpha=.4)
    paint(ax,geoms[(20,0)].intersection(region),'#117849','#117849')
    for x in [-10520,-10470,-9720,-9620]:ax.axvline(x,color='#687681',ls='--',lw=.7)
    for left,right,label in [(-10520,-10470,'50 µm接续'),(-10470,-9720,'MUX 本体：750 µm'),(-9720,-9620,'100 µm接续')]:
        ax.annotate('',xy=(right,1665),xytext=(left,1665),arrowprops={'arrowstyle':'<->','color':'#25475a'})
        ax.text((left+right)/2,1668,label,ha='center',fontsize=10 if right-left>100 else 8)
    ax.text(-10100,1635,'上方 MUX',ha='center',color='#15557e',fontsize=11)
    ax.text(-10100,1599,'下方 MUX',ha='center',color='#15557e',fontsize=11)
    ax.annotate('你标出的窗口接头',xy=(-10470,1586),xytext=(-10740,1617),
                arrowprops={'arrowstyle':'->','color':'#be402c'},color='#be402c',fontsize=10)
    ax.text(-10200,1620,'淡红色：GSG金属；绿色：LN1核心',fontsize=9,color='#784444')
    ax.set(xlim=(-10780,-9490),ylim=(1570,1683),xlabel='全局 x（µm）',ylabel='全局 y（µm）',
           title='左侧两只 MUX 的准确位置（实际版图几何；纵横不等比例）')
    ax.grid(alpha=.12)
    sections={}
    for j,(key,label) in enumerate([((20,1),'20/1：LN1刻蚀绘制窗口'),((10,2),'10/2：SiN挖除窗口')]):
        ax=fig.add_subplot(grid[1,j]);region=box(-10538,1574,-10410,1604)
        paint(ax,geoms[key].intersection(region),'#d9e7f5','#507797',origin=(-10470,1590))
        paint(ax,geoms[(20,0)].intersection(region),'#148555','#148555',origin=(-10470,1590))
        ax.axvline(0,color='#ba4932',ls='--',lw=1)
        section=[]
        for x in [-10470.01,-10469.99]:
            cut=geoms[key].intersection(LineString([(x,1574),(x,1604)]))
            section.append({'x_um':x,'y_min_um':cut.bounds[1],'y_max_um':cut.bounds[3],'width_um':cut.length})
        step=abs(section[0]['y_min_um']-section[1]['y_min_um']);yb=(section[0]['y_min_um']+section[1]['y_min_um'])/2-1590
        ax.annotate(f'边界台阶 {step:.2f} µm',xy=(0,yb),xytext=(12,-12.8),
                    arrowprops={'arrowstyle':'->','color':'#ba4932'},color='#ba4932',fontsize=10)
        ax.text(-35,10.8,f"接续侧：{section[0]['width_um']:.2f} µm",ha='center',fontsize=9)
        ax.text(31,10.8,f"本体侧：{section[1]['width_um']:.2f} µm",ha='center',fontsize=9)
        ax.set(xlim=(-68,60),ylim=(-16,14),xlabel='相对接头 x（µm）；0对应全局−10470 µm',
               ylabel='相对主路 y（µm）',title=label)
        ax.grid(alpha=.12);sections[str(key)]={'sections':section,'edge_step_um':step}
    output=a.output_dir/'MUX位置与窗口台阶.png'
    fig.savefig(output,dpi=170);plt.close(fig)
    report={'source_geometry':str(a.directory/'合并回读_四程几何.json'),
            'left_mux_body_x_um':[-10470,-9720],'left_external_x_um':[-10520,-10470],
            'left_active_connection_x_um':[-9720,-9620],'lower_window_steps':sections,
            'GDS_modified':False,'interpretation':'窗口阶跃与核心断线是不同检查；本图不把20/1当LN2材料边界'}
    with (a.output_dir/'窗口台阶核对.json').open('x',encoding='utf-8') as f:json.dump(report,f,ensure_ascii=False,indent=2)
    print(output.resolve())


if __name__=='__main__':main()
