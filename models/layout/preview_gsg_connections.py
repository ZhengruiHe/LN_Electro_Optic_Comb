"""显示新GSG工作稿，并从实际GDS导出的多边形绘制带中文标注的接口图。"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PatchPolygon
from shapely.geometry import Polygon,box
from klink import KLinkClient


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--show',action='store_true')
    a=parser.parse_args();directory=a.directory.resolve()
    if (directory/'GSG连接与版图检查.json').is_file():
        report=json.loads((directory/'GSG连接与版图检查.json').read_text(encoding='utf-8'))
    else:
        repaired=json.loads((directory/'分叉窗口与GSG检查.json').read_text(encoding='utf-8'))
        report={**repaired['source_report'],'output_gds':repaired['output_gds'],'top_cell':repaired['top_cell']}
    geometry=json.loads((directory/'版图预览几何.json').read_text(encoding='utf-8'))
    pitch=report['pitch_um']
    views=[('左端GSG',[1120,-255,2240,265],'RF_IN'),('右端GSG',[16760,-255,17880,265],'RF_OUT')]
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    colors={(10,2):'#d5ecf7',(20,1):'#b3dfca',(20,0):'#aa1748',(42,0):'#dcb451'}
    zorders={(10,2):1,(20,1):2,(42,0):3,(20,0):4}
    for title,bounds,interface in views:
        output=directory/(title+'_连接说明.png')
        if output.exists():raise FileExistsError(output)
        fig,ax=plt.subplots(figsize=(12,5.8),layout='constrained')
        window=box(*bounds)
        for layer in geometry:
            key=tuple(layer['layer'])
            for record in layer['polygons']:
                p=Polygon(record['exterior'],record['holes']).intersection(window)
                if p.is_empty:continue
                geoms=[p] if p.geom_type=='Polygon' else [g for g in p.geoms if g.geom_type=='Polygon']
                for geom in geoms:
                    ax.add_patch(PatchPolygon(list(geom.exterior.coords),facecolor=colors[key],
                                             edgecolor=colors[key],linewidth=.6,zorder=zorders[key]))
                    for hole in geom.interiors:
                        ax.add_patch(PatchPolygon(list(hole.coords),facecolor='white',edgecolor=colors[key],linewidth=.4,zorder=zorders[key]))
        points=report['probe_contact_centers_um'][interface]
        for net,(x,y) in points.items():
            ax.plot(x,y,'o',color='#1d3354',markersize=7,zorder=6)
            ax.text(x+42,y+12,'S 信号' if net=='S' else 'G 地',color='#1d3354',fontsize=11,zorder=6)
        px=points['S'][0];direction=1 if interface=='RF_IN' else -1
        ax.annotate('外接GSG探针',xy=(px,0),xytext=(px-direction*180,0),
                    ha='right' if direction==1 else 'left',va='center',fontsize=11,
                    arrowprops={'arrowstyle':'->','color':'#1d3354'},color='#1d3354',zorder=8)
        dx=px-direction*75
        for y0,y1 in [(0,pitch),(-pitch,0)]:
            ax.annotate('',xy=(dx,y0),xytext=(dx,y1),arrowprops={'arrowstyle':'<->','color':'#1d3354'},zorder=8)
            ax.text(dx-direction*5,(y0+y1)/2,f'{pitch:g} µm',rotation=90,
                    va='center',ha='right' if direction==1 else 'left',fontsize=10,zorder=8)
        ax.text((bounds[0]+bounds[2])/2,-232,'金色：M1金属  ·  红色：LN1光波导  ·  圆点：针尖接触中心',
                ha='center',fontsize=10)
        ax.set(xlim=(bounds[0],bounds[2]),ylim=(bounds[1],bounds[3]),aspect='equal',
               xlabel='版图x (µm)',ylabel='版图y (µm)',
               title=f'{title}：{pitch:g} µm针距，50 µm信号焊盘接43 µm主干')
        ax.grid(alpha=.15)
        fig.suptitle('几何连接已检查；阻抗、回波和实际探针头避让尚待验证',fontsize=10,color='#555')
        fig.savefig(output,dpi=170);plt.close(fig)
    if a.show:
        c=KLinkClient();c.connect()
        try:
            current=c.layout_info('full')
            if not current.get('file') or Path(current['file']).resolve()!=Path(report['output_gds']).resolve():
                c.layout_show_file(report['output_gds'],mode='new')
            c.call('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})
            for layer,color in [('10/2','#b9deee'),('20/1','#76c4a0'),('20/0','#b71c46'),('42/0','#d6ab43')]:
                c.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,
                                         'dither_pattern':0 if layer in ['20/0','42/0'] else 1,'line_width':1})
            c.call('layer.set_visible',{'layers':['101/0','56/30','100/30'],'visible':False})
            c.call('view.hier_levels',{'min':0,'max':10})
            for title,bounds,_ in views:
                c.call('view.screenshot',{'mode':'path','path':str(directory/(title+'_实际GDS.png')),
                       'bbox_um':bounds,'width_px':1700,'height_px':800})
            c.call('view.zoom_fit',{})
            c.call('view.screenshot',{'mode':'path','path':str(directory/'双端GSG_整体版图.png'),
                    'width_px':2400,'height_px':700})
            c.call('view.zoom_box',{'bbox_um':views[0][1]})
        finally:c.close()
    print(directory)


if __name__=='__main__':main()
