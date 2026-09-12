"""绘制已合并GDS的回读图形，不修改GDS。"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from search_fixed_block_placement import draw_layer


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--overwrite',action='store_true',help='仅重绘本脚本生成的预览PNG，不修改GDS')
    a=p.parse_args();d=a.directory
    old=json.loads(Path('results/layout/_历史归档/迭代版本_20260912/合并YSJ_版图确认_20260909_v1/YSJ源版图几何.json').read_text(encoding='utf-8'))
    new=json.loads((d/'合并回读_四程几何.json').read_text(encoding='utf-8'))
    manifest=json.loads((d/'合并清单.json').read_text(encoding='utf-8'))
    centers=sorted(tuple(c['center_um']) for c in manifest['crossings'])
    horizontal=manifest.get('layout_variant')=='horizontal_access'
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    views=[('合并版图_全图.png',[-11200,-2050,11200,2050],(18,4.5),'15 mm四程 + YSJ：同一21.8 × 3.8 mm block（功能版图候选）'),
           ('左侧延时与输入输出.png',[-10960,-1790,-6800,1910],(11,9),'左侧延时与输入输出：C1–C7为PDK交叉器'),
           ('右侧MUX与欧拉折返.png',[5200,1420,6800,1910],(13,5),'右侧MUX及四条欧拉折返：未改变YSJ结构')]
    if horizontal:
        views[1]=('左侧延时与输入输出.png',[-10960,-1790,-2500,1910],(15,7),f'横向接入与独立周期回路：C1–C{len(centers)}为PDK交叉器')
        views.append(('左侧横向接入细节.png',[-10930,980,-9000,1700],(14,6),'横向S接入与分离的欧拉折返（实际GDS轮廓）'))
    for filename,bounds,size,title in views:
        output=d/filename
        if output.exists() and not a.overwrite:raise FileExistsError(output)
        fig,ax=plt.subplots(figsize=size,layout='constrained')
        fig.get_layout_engine().set(rect=(0,.06,1,.94))
        for key,old_color,new_color in [((100,0),'#999999','#c87b1c'),((20,0),'#397758','#be2355'),((42,0),'#ac852f','#2972b8')]:
            draw_layer(ax,old,key,0,0,old_color)
            draw_layer(ax,new,key,0,0,new_color)
        ax.add_patch(Rectangle((-10.9,-1.9),21.8,3.8,fill=False,edgecolor='#555',linestyle='--'))
        if filename=='左侧延时与输入输出.png':
            for i,(x,y) in enumerate(centers,1):
                ax.text((x+40)/1000,(y+15)/1000,f'C{i}',fontsize=8,color='#a05b00')
            ports=manifest.get('edge_coupler_internal_ports_um',[[-10380.,-1700.],[-10380.,-500.]])
            for label,y,incoming in [('光输入',ports[0][1],True),('光输出',ports[1][1],False)]:
                start,end=((-10.87,-10.63) if incoming else (-10.63,-10.87))
                ax.annotate('',xy=(end,y/1000),xytext=(start,y/1000),
                            arrowprops={'arrowstyle':'->','color':'#555'})
                ax.text(-10.53,(y-65)/1000,label,fontsize=10,color='#555')
        ax.set(xlim=(bounds[0]/1000,bounds[2]/1000),ylim=(bounds[1]/1000,bounds[3]/1000),
               aspect='equal',xlabel='x（mm）',ylabel='y（mm）',title=title)
        ax.grid(alpha=.12)
        fig.text(.5,.012,'红/蓝：四程波导与电极；绿/棕：YSJ原结构。几何候选不等于完整光电性能或流片签核。',ha='center',fontsize=9,color='#555')
        fig.savefig(output,dpi=180);plt.close(fig)
    print(d.resolve())


if __name__=='__main__':main()
