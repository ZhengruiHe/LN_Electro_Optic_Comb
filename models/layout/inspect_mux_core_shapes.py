"""只读列出MUX直接20/0图形，辅助固定宽度窗口重建。"""
import json,pya
l=pya.Layout();l.read(gds_path);c=l.cell(cell_name)
rows=[]
for i,s in enumerate(c.shapes(l.layer(20,0)).each()):
    rows.append({'index':i,'type':'text' if s.is_text() else 'polygon' if s.is_polygon() else 'box' if s.is_box() else 'path' if s.is_path() else 'other',
                 'bbox_um':[s.dbbox().left,s.dbbox().bottom,s.dbbox().right,s.dbbox().top],
                 'text':s.text.string if s.is_text() else None})
with open(report_path,'x',encoding='utf-8') as f:json.dump(rows,f,ensure_ascii=False,indent=2)
