"""导出可提交的小型自有快照；不含PDK/黑盒或完整私有合并版图，不运行电磁求解。"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'results/layout/_历史归档/迭代版本_20260912/迭代版本归档清单_20260912.json'
GDS=ROOT/'results/layout/YSJ合并15mm_M04_KLink_v7_HTR有效电阻修正/YSJ与四程15mm_M04_O型输入_HTR有效50R_候选.gds'
CELLS=['MUX_TE01_A70_M04_A100_E50_BODY750_DRAFT','GSG150_LAUNCH_O_STYLE_50UM_DRAFT',
       'RF_TERMINATION_50R_RS20_EST_DRAFT','T_GSG_15MM_RECTANGULAR_BASELINE',
       'DUAL_RAIL_ACTIVE_LN_PARALLEL_WINDOWS_UNVERIFIED']
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def sanitized(x):
    if isinstance(x,dict): return {k:sanitized(v) for k,v in x.items()}
    if isinstance(x,list): return [sanitized(v) for v in x]
    if isinstance(x,str):
        x=x.replace(str(ROOT),'[PROJECT_ROOT]').replace(ROOT.as_posix(),'[PROJECT_ROOT]')
        if ':\\' in x or ':/Program Files/' in x: return '[外部本地文件，未包含]'
    return x
EXPORT=r"""
import json, pya
src=pya.Layout(); src.read(payload['source']); rows=[]
for name in payload['cells']:
    c=src.cell(name)
    if c is None or list(c.each_inst()): raise RuntimeError('只允许明确列出的无子实例自定义单元')
    prohibited={(70,30),(56,30),(100,0),(100,30)}
    if any(not c.shapes(li).is_empty() and (src.get_info(li).layer,src.get_info(li).datatype) in prohibited
           for li in src.layer_indices()):
        raise RuntimeError('单元含黑盒/PDK标记层，拒绝导出')
    dst=pya.Layout(); dst.dbu=src.dbu
    target=dst.create_cell(name); target.copy_tree(c)
    if len(list(dst.each_cell()))!=1: raise RuntimeError('意外复制了子树')
    path=payload['out']+'/'+name+'.gds'; dst.write(path)
    check=pya.Layout(); check.read(path)
    for li in src.layer_indices():
        if c.shapes(li).is_empty(): continue
        info=src.get_info(li); idx=check.find_layer(info.layer,info.datatype)
        if idx is None: raise RuntimeError('导出丢失图层')
        if not (pya.Region(c.shapes(li))^pya.Region(check.top_cell().shapes(idx))).is_empty():
            raise RuntimeError('导出几何不一致')
    rows.append({'cell':name,'path':path,'roundtrip_geometry_equal':True})
mux=src.cell(payload['cells'][0])
core=pya.Region(mux.shapes(src.layer(20,0))).merged()
clad=pya.Region(mux.shapes(src.layer(20,1))).merged()
etch=(clad-core).merged()
rules={'core_width_030':core.width_check(300).size(),'core_space_030':core.space_check(300).size(),
       'etch_width_030':etch.width_check(300).size(),'etch_space_030':etch.space_check(300).size()}
json.dumps({'cells':rows,'mux_local_markers':rules})
"""
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args(); out=a.output_dir.resolve()
    if not out.is_relative_to(ROOT) or out.exists(): raise ValueError('需要项目内不存在的新输出目录')
    manifest=json.loads(ARCHIVE.read_text(encoding='utf-8'))
    for row in manifest['moved']:
        d=Path(row['archive_path']).resolve()
        if not d.is_relative_to(ROOT/'results/layout/_历史归档'): raise ValueError('归档越界')
        for f in row['files']:
            path=(d/f['relative_path']).resolve()
            if not path.is_relative_to(d) or path.stat().st_size!=f['size'] or digest(path)!=f['sha256']:
                raise RuntimeError('归档文件校验失败：'+str(path))
    out.mkdir(parents=True); (out/'self_cells').mkdir()
    archive_rows=[{'source':Path(r['source_path']).relative_to(ROOT).as_posix(),
                   'archive':Path(r['archive_path']).relative_to(ROOT).as_posix(),
                   'files':r['files']} for r in manifest['moved']]
    def write(name,data): (out/name).write_text(json.dumps(sanitized(data),ensure_ascii=False,indent=2),encoding='utf-8',newline='\n')
    write('archive_index.json',{'directories':archive_rows,'files_hash_verified':True,'deleted_files':0})
    from klink import KLinkClient
    c=KLinkClient(); c.connect()
    try:
        before=digest(GDS)
        r=c.call('exec.python',{'code':'payload='+repr({'source':str(GDS),'out':str(out/'self_cells'),'cells':CELLS})+'\n'+EXPORT})
        if r.get('exception'): raise RuntimeError(r['exception'])
        exported=json.loads(r['return_value'])
        if digest(GDS)!=before: raise RuntimeError('源GDS变化')
    finally: c.close()
    write('mux_local_geometry_check.json',{'source_gds_sha256':before,'markers':exported['mux_local_markers'],
          'scope':'仅自定义MUX单元；非完整DRC','foundry_signoff':False})
    sources={
        'resistor_readback.json':GDS.parent/'电阻有效长度与接触回读检查.json',
        'Al_absorption_simplified.json':ROOT/'results/optical/金属吸收_Al局部_20260912_v1/金属吸收损耗结果.json',
        'TiN_absorption_simplified.json':ROOT/'results/optical/金属吸收_TiN局部_20260912_v1/金属吸收损耗结果.json'}
    for name,src in sources.items():
        data=json.loads(src.read_text(encoding='utf-8'))
        write(name,{'source_sha256':digest(src),'original_result':data,
                    'review_note':('简化截面筛查：刻蚀区SiO2填充有遗漏，非实际桥结构总损耗；用户接受用于当前阶段，非收敛签核。'
                                   if 'absorption' in name else '几何电阻按20 ohm/square用户假设，不是RF实测。')})
    touch=ROOT/'results/hfss/O_input_screening_v5/O_input_GSG150_8_12GHz.s4p'
    text=touch.read_text(encoding='utf-8')
    text='\n'.join(line.rstrip() for line in text.splitlines() if not line.startswith('!        File:'))+'\n'
    (out/'O_input_GSG150_8_12GHz.s4p').write_text(text,encoding='utf-8',newline='\n')
    nums=[list(map(float,l.split())) for l in text.splitlines() if l.strip() and not l.lstrip().startswith(('!','#'))]
    if '# GHz S MA' not in text or '! Port[2] = P1_CPW' not in text or '! Port[4] = P2_CPW' not in text:
        raise ValueError('端口顺序或Touchstone格式变化')
    records=[]
    for i in range(0,len(nums),4):
        f=nums[i][0]; matrix=[nums[i][1:]]+nums[i+1:i+4]
        records.append({'frequency_ghz':f,'S11_CPW_db':20*math.log10(matrix[1][2]),
                        'S21_CPW_db':20*math.log10(matrix[3][2])})
    write('hfss_input_cpw_summary.json',{'source_sha256':digest(touch),'records':records,
        'physical_ports':['P1_CPW','P2_CPW'],'scope':'321um输入+500um均匀CPW；不含15mm T电极和HTR终端',
        'rf_signoff':False,'correction':'早期-25.78/-0.085dB为Odd模，不用作GSG物理CPW结论。'})
    files=[{'path':f.relative_to(out).as_posix(),'bytes':f.stat().st_size,'sha256':digest(f)}
           for f in sorted(out.rglob('*')) if f.is_file()]
    write('snapshot_manifest.json',{'source_layout_sha256':before,'files':files,'pdk_included':False,
          'contains_full_chip':False,'scope':'自定义单元小库、轻量数值结果与归档索引；非生产签核包'})
    print(json.dumps({'files':len(files)+1,'bytes':sum(x['bytes'] for x in files),'output':str(out)},ensure_ascii=False))
if __name__=='__main__': main()
