"""实际GDS回读窗口的局部几何预览，不修改GDS。"""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as Patch
from shapely.geometry import box
from search_fixed_block_placement import geometry
ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
d=json.loads((a.directory/'四程实际掩膜预览.json').read_text(encoding='utf-8'))
fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')
for ax,(x,y,key) in zip(axes,[(-10676.,1876.7,(20,1)),(-10682.8,1797.7,(10,2))]):
    region=box(x-8,y-8,x+8,y+8)
    for k,c in [(key,'#bcd8ef'),((20,0),'#226c32')]:
        g=geometry(d,k).intersection(region)
        for p in g.geoms if hasattr(g,'geoms') else [g]:
            if p.is_empty or p.geom_type!='Polygon':continue
            ax.add_patch(Patch(p.exterior.coords,facecolor=c,edgecolor='#345',lw=.5))
            for h in p.interiors:ax.add_patch(Patch(h.coords,facecolor='white',edgecolor='#345',lw=.5))
    ax.set(xlim=(x-8,x+8),ylim=(y-8,y+8),aspect='equal',title=str(key));ax.grid(alpha=.2)
fig.savefig(a.output,dpi=160)
