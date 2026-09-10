"""在KLink中新开版图页并保存实际GDS截图；另以原生几何画拉锥尺寸核对图。"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from klink import KLinkClient
from shapely.geometry import Polygon, box
from shapely.ops import unary_union


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    args = parser.parse_args()
    report = json.loads((args.directory/'版图检查.json').read_text(encoding='utf-8'))
    client = KLinkClient()
    client.connect()
    try:
        for method,params in [('layout.show_file',{'path':report['working_gds'],'mode':'new'}),
                              ('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})]:
            client.call(method,params)
        for layer,color in (('10/2','#C6E8FB'),('20/1','#71BBA9'),('20/0','#BE2654'),
                            ('42/0','#D5A12C'),('100/0','#766A85'),('70/30','#714191')):
            client.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,
                                          'dither_pattern':0 if layer in ('20/0','42/0','70/30') else 1,'line_width':1})
        client.call('layer.set_visible',{'layers':['101/0','100/30','56/30'],'visible':False})
        client.call('view.hier_levels',{'min':0,'max':10})
        for name,bbox,w,h in [('整体版图.png',[-620,-560,19320,1040],2400,500),
                              ('黑盒交叉与四端拉锥.png',[160,110,440,390],1400,1400),
                              ('输入黑盒与交叉接续.png',[-540,80,560,440],2200,800)]:
            client.call('view.screenshot',{'mode':'path','path':str((args.directory/name).resolve()),
                                           'bbox_um':bbox,'width_px':w,'height_px':h})
        client.call('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})
    finally:
        client.close()

    preview = json.loads((args.directory/'预览几何.json').read_text(encoding='utf-8'))
    layer = next(v for v in preview if v['layer']==[20,0])
    core = unary_union([Polygon(v['exterior'],v['holes']) for v in layer['polygons']])
    length = report['taper_length_um']
    start,end = 277.5-length,277.5
    local = core.intersection(box(start,248,end,252))
    if local.geom_type != 'Polygon':
        raise ValueError('实际GDS西侧拉锥不是单连通区域')
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','DejaVu Sans'],'axes.unicode_minus':False})
    fig,ax = plt.subplots(figsize=(11,3.8),layout='constrained')
    pts = list(local.exterior.coords)
    ax.fill([p[0]-start for p in pts],[p[1]-250 for p in pts],color='#BE2654',alpha=.85)
    ax.axhline(0,color='#666666',lw=.6,ls='--')
    ax.set(xlim=(-5,length+5),ylim=(-1.15,1.15),xlabel='沿拉锥方向的距离 (µm)',ylabel='横向 (µm)',
           title='实际 GDS 西侧 LN1 拉锥：0.7 → 1.2 µm，长度 '+str(length)+' µm')
    for x,width,label in ((0,.7,'接原路由'),(length,1.2,'接黑盒端面')):
        ax.annotate('',xy=(x,width/2),xytext=(x,-width/2),arrowprops={'arrowstyle':'<->','color':'#333333'})
        ax.text(x,.82,f'{label}\n{width:g} µm',ha='left' if x==0 else 'right',va='center')
    fig.suptitle('仅展示绘制轮廓，非光学验证；横纵比例不同',fontsize=10,color='#555555')
    fig.savefig(args.directory/'拉锥实际GDS尺寸核对.png',dpi=180)
    plt.close(fig)
    print(args.directory.resolve())


if __name__ == '__main__':
    main()
