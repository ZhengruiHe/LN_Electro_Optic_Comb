"""只读检查一份失败/候选GDS中各单元的窗口间距标记。"""
import json
import pya
with open(payload_file,encoding='utf-8') as f:p=json.load(f)
l=pya.Layout();l.read(gds_path if 'gds_path' in globals() else p['ours_gds'])
rows=[]
for cell in l.each_cell():
    if not cell.name.startswith(('MUX','EO4P')):continue
    core=pya.Region(cell.begin_shapes_rec(l.layer(20,0))).merged()
    clad=pya.Region(cell.begin_shapes_rec(l.layer(20,1))).merged()
    etch=clad-core
    marks=etch.space_check(300);width=etch.width_check(300)
    trench=pya.Region(cell.begin_shapes_rec(l.layer(10,2))).merged();twidth=trench.width_check(200)
    rows.append({'cell':cell.name,'core_area_um2':core.area()*l.dbu*l.dbu,
                 'clad_area_um2':clad.area()*l.dbu*l.dbu,'markers':marks.size(),
                 'boxes_um':[[e.bbox().left*l.dbu,e.bbox().bottom*l.dbu,e.bbox().right*l.dbu,e.bbox().top*l.dbu] for e in marks.each()],
                 'width_markers':width.size(),'width_boxes_um':[[e.bbox().left*l.dbu,e.bbox().bottom*l.dbu,e.bbox().right*l.dbu,e.bbox().top*l.dbu] for e in width.each()],
                 'trench_width_markers':twidth.size(),'trench_width_boxes_um':[[e.bbox().left*l.dbu,e.bbox().bottom*l.dbu,e.bbox().right*l.dbu,e.bbox().top*l.dbu] for e in twidth.each()]})
print(json.dumps(rows,ensure_ascii=False,indent=2))
if 'report_path' in globals():
    with open(report_path,'x',encoding='utf-8') as f:json.dump(rows,f,ensure_ascii=False,indent=2)
