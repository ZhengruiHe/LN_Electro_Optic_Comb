"""四程版图的部分欧拉弯中心线：纯几何，不写GDS、不声称弯曲损耗通过。"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

import numpy as np

from generate_four_pass_gds import polyline_length

EULER_FRACTION = 0.40
MIN_RADIUS_UM = 80.0


@dataclass
class BendInfo:
    angle_deg: float
    minimum_radius_um: float
    euler_fraction: float
    centerline_length_um: float
    displacement_x_um: float
    displacement_y_um: float
    start_heading_deg: float
    end_heading_deg: float


def partial_euler_local(delta_angle: float, minimum_radius: float,
                        euler_fraction: float = EULER_FRACTION,
                        maximum_step_um: float = 0.5) -> tuple[np.ndarray, dict]:
    """从原点沿+x出发，返回曲率连续的部分欧拉弯折线。"""
    angle = abs(delta_angle)
    if not (0 < angle <= math.pi) or minimum_radius < MIN_RADIUS_UM:
        raise ValueError("弯曲角必须在0到180度，最小曲率半径不得低于80um")
    if not 0 < euler_fraction <= 1 or maximum_step_um <= 0:
        raise ValueError("欧拉比例或采样步长无效")
    sign = 1.0 if delta_angle > 0 else -1.0
    ramp_length = minimum_radius * euler_fraction * angle
    circle_length = minimum_radius * (1.0-euler_fraction) * angle
    sections = (("ramp_up",ramp_length),("circle",circle_length),("ramp_down",ramp_length))
    x=y=heading=0.0
    points=[[0.0,0.0]]
    for name,length in sections:
        if length <= 0:
            continue
        steps=max(1,math.ceil(length/maximum_step_um))
        ds=length/steps
        for index in range(steps):
            s=(index+0.5)*ds
            if name=="ramp_up":
                curvature=s/(minimum_radius*ramp_length)
            elif name=="circle":
                curvature=1.0/minimum_radius
            else:
                curvature=(1.0-s/ramp_length)/minimum_radius
            curvature*=sign
            middle=heading+0.5*curvature*ds
            x+=math.cos(middle)*ds
            y+=math.sin(middle)*ds
            heading+=curvature*ds
            points.append([x,y])
    expected=sign*angle
    if not math.isclose(heading,expected,abs_tol=2e-10):
        raise RuntimeError("欧拉弯积分角度不闭合")
    return np.asarray(points), {
        "angle_deg":math.degrees(delta_angle),"minimum_radius_um":minimum_radius,
        "euler_fraction":euler_fraction,"centerline_length_um":2*ramp_length+circle_length,
        "displacement_x_um":x,"displacement_y_um":y,"end_heading_rad":heading,
        "maximum_curvature_per_um":1.0/minimum_radius,
    }


def append_euler(points: list[list[float]], heading: float, delta_angle: float,
                 minimum_radius: float, bend_log: list[dict], label: str,
                 euler_fraction: float = EULER_FRACTION) -> float:
    local,info=partial_euler_local(delta_angle,minimum_radius,euler_fraction=euler_fraction)
    start=np.asarray(points[-1],dtype=float)
    rotation=np.asarray([[math.cos(heading),-math.sin(heading)],
                         [math.sin(heading), math.cos(heading)]])
    transformed=local[1:]@rotation.T+start
    points.extend(transformed.tolist())
    end_heading=heading+delta_angle
    row=BendInfo(info["angle_deg"],minimum_radius,euler_fraction,
                 info["centerline_length_um"],float(transformed[-1,0]-start[0]),
                 float(transformed[-1,1]-start[1]),math.degrees(heading),
                 math.degrees(end_heading))
    bend_log.append({"label":label,**asdict(row)})
    return end_heading


def bend_displacement(angle: float, radius: float) -> tuple[float,float]:
    # 必须与append_euler使用相同积分步长，否则理论落点与实际落点会出现
    # 亚纳米差，随后被长直段放大成数学上的斜直线。
    _,info=partial_euler_local(angle,radius)
    return info["displacement_x_um"],info["displacement_y_um"]


def radius_for_semicircle_displacement(target_displacement_um: float,
                                       minimum_radius_um: float = MIN_RADIUS_UM) -> float:
    """反解默认离散欧拉半圆的纵向位移，保证端口高度精确闭合。"""
    if target_displacement_um<=0:
        raise ValueError("目标半圆位移必须为正")
    low=minimum_radius_um
    high=max(low*2.0,target_displacement_um)
    while bend_displacement(math.pi,high)[1]<target_displacement_um:
        high*=2.0
    if bend_displacement(math.pi,low)[1]>target_displacement_um+1e-12:
        raise ValueError("目标位移要求的曲率半径低于80um")
    for _ in range(60):
        middle=.5*(low+high)
        if bend_displacement(math.pi,middle)[1]<target_displacement_um:
            low=middle
        else:
            high=middle
    return .5*(low+high)


def normalized_bend_displacement(angle: float) -> tuple[float,float]:
    """用满足PDK的80um弯计算线性归一化位移，避免构造非物理1um弯。"""
    x,y=bend_displacement(angle,MIN_RADIUS_UM)
    return x/MIN_RADIUS_UM,y/MIN_RADIUS_UM


def euler_return_route(*, x_right: float, y_start: float, x_left: float,
                       y_end: float, target_length: float, outward_sign: int,
                       start_radius: float, middle_radius: float = MIN_RADIUS_UM,
                       launch_extension: float = 5.0,
                       end_extension: float = 5.0) -> tuple[list[list[float]],dict]:
    if outward_sign not in (-1,1):
        raise ValueError("outward_sign只能取-1或1")
    sign=float(outward_sign)
    _,d_start=bend_displacement(math.pi,start_radius)
    _,d_middle=bend_displacement(math.pi,middle_radius)
    lane1=y_start+sign*d_start
    lane2=lane1+sign*d_middle
    lane3=lane2+sign*d_middle
    final_radius=abs(lane3-y_end)/normalized_bend_displacement(math.pi)[1]
    if min(start_radius,middle_radius,final_radius)<MIN_RADIUS_UM-1e-8:
        raise ValueError("欧拉回路出现小于80um的最小曲率半径")
    arc_length=sum(r*math.pi*(1+EULER_FRACTION)
                   for r in (start_radius,middle_radius,middle_radius,final_radius))
    direct=x_right-x_left
    backtrack=.5*(target_length-direct-2*launch_extension-2*end_extension-arc_length)
    if backtrack<=0:
        raise ValueError("目标回路长度不足以容纳欧拉弯")
    x_turn=x_left+max(1300.0,1.5*final_radius+300.0)
    x_second=x_turn+backtrack
    if x_second+1.5*middle_radius>=x_right+launch_extension-1.5*start_radius:
        raise ValueError("欧拉蛇形横向空间不足")
    points=[[x_right,y_start]]
    if launch_extension:
        points.append([x_right+launch_extension,y_start])
    bends=[]
    heading=0.0
    heading=append_euler(points,heading,sign*math.pi,start_radius,bends,"右侧折返1")
    points.append([x_turn,lane1])
    heading=append_euler(points,heading,-sign*math.pi,middle_radius,bends,"左侧折返2")
    points.append([x_second,lane2])
    heading=append_euler(points,heading,sign*math.pi,middle_radius,bends,"右侧折返3")
    final_x=x_left-end_extension
    points.append([final_x,lane3])
    heading=append_euler(points,heading,sign*math.pi,final_radius,bends,"左侧落回4")
    if end_extension:
        points.append([x_left,y_end])
    actual=polyline_length(points)
    if not math.isclose(actual,target_length,abs_tol=2e-3):
        raise RuntimeError(f"欧拉回路长度回标失败：{actual} vs {target_length}")
    return points,{"target_length_um":target_length,"actual_centerline_length_um":actual,
                   "euler_fraction":EULER_FRACTION,"minimum_radius_um":min(b["minimum_radius_um"] for b in bends),
                   "start_radius_um":start_radius,"middle_radius_um":middle_radius,"end_radius_um":final_radius,
                   "horizontal_backtrack_um":backtrack,"outer_lane_y_um":lane3,"bends":bends}


def euler_compact_finish_return_route(*, x_right: float, y_start: float,
                                      x_left: float, y_end: float,
                                      target_length: float, outward_sign: int,
                                      middle_radius: float = MIN_RADIUS_UM,
                                      launch_extension: float = 5.0,
                                      end_extension: float = 5.0,
                                      finish_excursion_um: float = 655.0
                                      ) -> tuple[list[list[float]],dict]:
    """六个欧拉U弯的折返。

    最后三个U弯把外侧路线分三级落回端口，避免单个大半径终弯向左
    占用输入光路区域。该函数只处理中心线几何和长度回标。
    """
    if outward_sign not in (-1,1):
        raise ValueError("outward_sign只能取-1或1")
    sign=float(outward_sign)
    _,d_middle=bend_displacement(math.pi,middle_radius)
    target_start_displacement=d_middle+sign*(y_end-y_start)
    start_radius=radius_for_semicircle_displacement(
        target_start_displacement,middle_radius)
    if min(start_radius,middle_radius)<MIN_RADIUS_UM-1e-8:
        raise ValueError("紧凑落回结构出现小于80um的最小曲率半径")
    _,d_start=bend_displacement(math.pi,start_radius)
    lane1=y_start+sign*d_start
    lane2=lane1+sign*d_middle
    lane3=lane2+sign*d_middle
    lane4=lane3-sign*d_middle
    lane5=lane4-sign*d_middle
    calculated_end=lane5-sign*d_middle
    if not math.isclose(calculated_end,y_end,abs_tol=2e-3):
        raise RuntimeError("紧凑落回结构的纵向端点不闭合")

    arc_length=(start_radius+5*middle_radius)*math.pi*(1+EULER_FRACTION)
    direct=x_right-x_left
    backtrack=.5*(target_length-direct-2*launch_extension-2*end_extension
                  -2*finish_excursion_um-arc_length)
    if backtrack<=0:
        raise ValueError("目标回路长度不足以容纳六弯欧拉结构")
    final_x=x_left-end_extension
    finish_right_x=final_x+finish_excursion_um
    x_turn=x_left+1300.0
    x_second=x_turn+backtrack
    if x_second<=finish_right_x+2.0*middle_radius:
        raise ValueError("主蛇形与三级落回没有足够横向间隔")
    if x_second+1.5*middle_radius>=x_right+launch_extension-1.5*start_radius:
        raise ValueError("六弯欧拉蛇形横向空间不足")

    points=[[x_right,y_start]]
    if launch_extension:
        points.append([x_right+launch_extension,y_start])
    bends=[]; heading=0.0
    heading=append_euler(points,heading,sign*math.pi,start_radius,bends,"右侧折返1")
    lane1=points[-1][1]
    points.append([x_turn,lane1])
    heading=append_euler(points,heading,-sign*math.pi,middle_radius,bends,"左侧折返2")
    lane2=points[-1][1]
    points.append([x_second,lane2])
    heading=append_euler(points,heading,sign*math.pi,middle_radius,bends,"右侧折返3")
    lane3=points[-1][1]
    points.append([final_x,lane3])
    heading=append_euler(points,heading,sign*math.pi,middle_radius,bends,"左侧落回4")
    lane4=points[-1][1]
    points.append([finish_right_x,lane4])
    heading=append_euler(points,heading,-sign*math.pi,middle_radius,bends,"右侧落回5")
    lane5=points[-1][1]
    points.append([final_x,lane5])
    heading=append_euler(points,heading,sign*math.pi,middle_radius,bends,"左侧落回6")
    if not math.isclose(points[-1][1],y_end,abs_tol=1e-9):
        raise RuntimeError(f"六弯欧拉端口高度回标失败：{points[-1][1]} vs {y_end}")
    if end_extension:
        points.append([x_left,points[-1][1]])
    actual=polyline_length(points)
    if not math.isclose(actual,target_length,abs_tol=2e-3):
        raise RuntimeError(f"六弯欧拉回路长度回标失败：{actual} vs {target_length}")
    return points,{"target_length_um":target_length,"actual_centerline_length_um":actual,
                   "euler_fraction":EULER_FRACTION,"minimum_radius_um":min(b["minimum_radius_um"] for b in bends),
                   "start_radius_um":start_radius,"middle_radius_um":middle_radius,
                   "horizontal_backtrack_um":backtrack,"finish_excursion_um":finish_excursion_um,
                   "outer_lane_y_um":lane3,"bends":bends}


def euler_loop2(target_length: float) -> tuple[list[list[float]],dict]:
    xr,xl,ys,ye=18150.0,850.0,26.6,-30.15
    launch,r0,rm,rq=400.0,270.0,80.0,80.0
    xc,xt=300.0,2150.0
    _,d0=bend_displacement(math.pi,r0)
    _,dm=bend_displacement(math.pi,rm)
    dq,_=bend_displacement(math.pi/2,rq)
    lane1=ys+d0; lane2=lane1+dm; lane3=lane2+dm
    vertical=lane3-ye-2*dq
    if vertical<=0:
        raise ValueError("loop2欧拉90度弯没有竖直连接空间")
    arc_length=(r0+2*rm)*math.pi*(1+EULER_FRACTION)+2*rq*math.pi/2*(1+EULER_FRACTION)
    fixed=xr-xl+2*launch+2*(xl-xc-dq)+vertical+arc_length
    backtrack=.5*(target_length-fixed)
    if backtrack<=0:
        raise ValueError("loop2目标长度不足")
    x_second=xt+backtrack
    if x_second+1.5*rm>=xr+launch-1.5*r0:
        raise ValueError("loop2横向空间不足")
    points=[[xr,ys],[xr+launch,ys]]
    bends=[]; heading=0.0
    heading=append_euler(points,heading,math.pi,r0,bends,"回路2折返1")
    lane1=points[-1][1]
    points.append([xt,lane1])
    heading=append_euler(points,heading,-math.pi,rm,bends,"回路2折返2")
    lane2=points[-1][1]
    points.append([x_second,lane2])
    heading=append_euler(points,heading,math.pi,rm,bends,"回路2折返3")
    lane3=points[-1][1]
    vertical=lane3-ye-2*dq
    if vertical<=0:
        raise ValueError("loop2欧拉90度弯没有竖直连接空间")
    points.append([xc+dq,lane3])
    heading=append_euler(points,heading,math.pi/2,rq,bends,"回路2转竖直4")
    points.append([points[-1][0],ye+dq])
    heading=append_euler(points,heading,math.pi/2,rq,bends,"回路2转水平5")
    if not math.isclose(points[-1][1],ye,abs_tol=1e-9):
        raise RuntimeError("loop2端口高度回标失败")
    points.append([xl,points[-1][1]])
    actual=polyline_length(points)
    if not math.isclose(actual,target_length,abs_tol=2e-3):
        raise RuntimeError(f"loop2欧拉长度回标失败：{actual} vs {target_length}")
    return points,{"target_length_um":target_length,"actual_centerline_length_um":actual,
                   "euler_fraction":EULER_FRACTION,"minimum_radius_um":min(b["minimum_radius_um"] for b in bends),
                   "horizontal_backtrack_um":backtrack,"outer_lane_y_um":lane3,"bends":bends}


def euler_input_route(bend_start_x: float = 445.0,
                      euler_fraction: float = EULER_FRACTION) -> tuple[list[list[float]],dict]:
    y0,y1=250.0,30.15
    radius=80.0
    _,quarter=partial_euler_local(math.pi/2,radius,euler_fraction=euler_fraction)
    dq=quarter['displacement_x_um']
    vertical=y0-y1-2*dq
    if vertical<0:
        raise ValueError("输入端高度不足以容纳80um欧拉S弯")
    if not 422.5 <= bend_start_x <= 600:
        raise ValueError("输入欧拉S弯起点必须位于交叉器东端与MUX之间")
    points=[[0.0,y0],[bend_start_x,y0]]
    bends=[]; heading=0.0
    heading=append_euler(points,heading,-math.pi/2,radius,bends,"输入转下1",euler_fraction)
    points.append([points[-1][0],y1+dq])
    heading=append_euler(points,heading,math.pi/2,radius,bends,"输入转平2",euler_fraction)
    if not math.isclose(points[-1][1],y1,abs_tol=1e-9):
        raise RuntimeError("输入欧拉S弯端口高度回标失败")
    if points[-1][0] >= 850:
        raise ValueError("输入欧拉S弯没有留出接MUX的直段")
    points.append([850.0,points[-1][1]])
    return points,{"actual_centerline_length_um":polyline_length(points),"vertical_straight_um":vertical,
                   "euler_fraction":euler_fraction,"minimum_radius_um":radius,"bends":bends}


def build_euler_routes(target_lengths: list[float], input_bend_start_x: float = 445.0,
                       input_euler_fraction: float = EULER_FRACTION) -> tuple[dict,dict]:
    if len(target_lengths)!=3:
        raise ValueError("需要三个回路目标长度")
    loop1,info1=euler_compact_finish_return_route(
        x_right=18150,x_left=850,y_start=29.85,y_end=33.4,
        target_length=target_lengths[0],outward_sign=1)
    loop2,info2=euler_loop2(target_lengths[1])
    loop3,info3=euler_compact_finish_return_route(
        x_right=18150,x_left=850,y_start=-29.85,y_end=-33.4,
        target_length=target_lengths[2],outward_sign=-1)
    incoming,input_info=euler_input_route(input_bend_start_x,input_euler_fraction)
    routes={"loop1":loop1,"loop2":loop2,"loop3":loop3,"input":incoming,
            "output":[[18150,-26.6],[18700,-26.6]]}
    return routes,{"loop1":info1,"loop2":info2,"loop3":info3,"input":input_info,
                   "total_euler_bends":sum(len(v["bends"]) for v in (info1,info2,info3,input_info))}
