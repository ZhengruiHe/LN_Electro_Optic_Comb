"""在独立KLayout页打开指定工作稿，导出整体和交叉器局部预览。"""
import argparse
from pathlib import Path
from klink import KLinkClient


def main():
    p = argparse.ArgumentParser()
    p.add_argument("gds", type=Path)
    p.add_argument("--top-cell", default="EO4P_10G_15MM_CUSTOM_CROSSING")
    a = p.parse_args()
    c = KLinkClient()
    c.connect()
    try:
        c.call("layout.show_file", {"path":str(a.gds.resolve()),"mode":"new"})
        c.call("view.show_cell", {"cell":a.top_cell,"zoom_fit":True})
        for layer, color in (("10/2","#C6E8FB"),("20/1","#71BBA9"),
                             ("20/0","#BE2654"),("42/0","#D5A12C")):
            c.call("layer.set_style", {"layer":layer,"fill_color":color,
                                         "frame_color":color,"dither_pattern":0,"line_width":1})
        c.call("layer.set_visible", {"layers":["101/0"],"visible":False})
        c.call("view.hier_levels", {"min":0,"max":10})
        for name, bbox, width, height in (
            ("整体版图.png", [-620,-560,19320,950],2400,600),
            ("左侧正交交叉与回路.png", [-50,-80,1000,980],1200,1250),
            ("自定义交叉器_版图局部.png", [275,225,325,275],1200,1200),
        ):
            path=a.gds.parent/name
            c.call("view.screenshot", {"mode":"path","path":str(path.resolve()),
                                        "bbox_um":bbox,"width_px":width,"height_px":height})
            print(path)
    finally:
        c.close()


if __name__ == "__main__":
    main()
