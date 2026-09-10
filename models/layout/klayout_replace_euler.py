"""由KLayout批处理调用：清除顶层旧无源路由并写入欧拉多边形清单。"""
import json
import pya


required=('source','output','manifest','source_top','output_top')
missing=[name for name in required if name not in globals()]
if missing:
    raise RuntimeError('缺少-rd参数：'+','.join(missing))

layout=pya.Layout()
layout.read(source)
cell=layout.cell(source_top)
if cell is None:
    raise RuntimeError('找不到源顶层单元：'+source_top)

data=json.loads(open(manifest,'r',encoding='utf-8').read())
for layer_name in ('10/2','20/0','20/1'):
    layer,datatype=(int(item) for item in layer_name.split('/'))
    layer_index=layout.layer(pya.LayerInfo(layer,datatype))
    cell.shapes(layer_index).clear()
    for polygon in data['polygons'][layer_name]:
        points=[pya.DPoint(float(x),float(y)) for x,y in polygon]
        cell.shapes(layer_index).insert(pya.DPolygon(points))

cell.name=output_top
layout.write(output)
