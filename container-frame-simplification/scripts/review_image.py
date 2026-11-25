"""Render grouped review images from the same CAD surfaces measured for G+C."""

from functools import lru_cache
import os
from pathlib import Path

import numpy as np
from PIL import Image, ImageColor, ImageDraw, ImageFont
from OCP.BRep import BRep_Tool
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopAbs import TopAbs_FACE
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS

from analysis_bundle import read_analysis
from common import load, save, write_csv
from intake import Model, indexed
from measure_geometry import gap
from raster import Camera, rasterize


ORANGE = "#e87929"
BLUE = "#3b7ec0"
GRAY = "#8e9caa"
INK = "#1c384d"
FRAME_ROLES = {"CAST","POST","POST_TWO","POST_ONE","TOP_RAIL","BOTTOM_RAIL","CROSS","DECK","FRAME_SEG","BOTTOM_LONG"}
PAGE_WIDTH = 1500
MARGIN = 26


@lru_cache(maxsize=24)
def font(size, bold=False):
    fonts = Path(os.environ.get("WINDIR","C:/Windows"))/"Fonts"
    return ImageFont.truetype(str(fonts/("msyhbd.ttc" if bold else "msyh.ttc")),size)


def label(draw, point, text, size=20, color=INK, bold=False, anchor=None):
    draw.text(point,text,font=font(size,bold),fill=color,anchor=anchor)


def triangles(shape, tolerance):
    BRepMesh_IncrementalMesh(shape,tolerance,False,.28,False)
    output = []
    for raw in indexed(shape,TopAbs_FACE):
        face = TopoDS.Face_s(raw)
        location = TopLoc_Location()
        mesh = BRep_Tool.Triangulation_s(face,location)
        if mesh is None:
            continue
        transform = location.Transformation()
        nodes = [mesh.Node(n).Transformed(transform).Coord() for n in range(1,mesh.NbNodes()+1)]
        output.extend([nodes[j-1] for j in mesh.Triangle(n).Get()] for n in range(1,mesh.NbTriangles()+1))
    return np.asarray(output,dtype=np.float32).reshape((-1,3,3))


def render(camera, solids, colors):
    verts = camera.project(np.concatenate(solids))*1.5
    sizes = [len(s) for s in solids]
    ids = np.repeat(np.arange(1,len(solids)+1,dtype=np.int32),sizes)
    palette = np.array([(0,0,0),*[ImageColor.getrgb(x) for x in colors]],dtype=np.uint8)
    return Image.fromarray(rasterize(np.ascontiguousarray(verts,dtype=np.float32),ids,palette,
                          round(camera.width*1.5),round(camera.height*1.5))).resize((camera.width,camera.height),Image.Resampling.LANCZOS)


def overview(base, frame, pending, examples, width, height):
    direction = [-.16,-.72,.66]
    bounds = np.asarray(base["frame_box"])
    guide = Camera(width,height,bounds.mean(axis=0),direction,1.,up=(0,0,1))
    relative = (np.concatenate(frame).reshape(-1,3)-guide.center) @ guide.axes.T
    low, high = relative.min(axis=0),relative.max(axis=0)
    guide.center += ((low+high)/2) @ guide.axes
    guide.scale = min((width-130)/(high[0]-low[0]),(height-90)/(high[1]-low[1]))
    view = render(guide,frame,[GRAY]*len(frame))
    draw = ImageDraw.Draw(view)
    for ident in pending:
        if "box" in base["parts"][ident]:
            xyz = np.asarray(base["parts"][ident]["box"]).mean(axis=0)
            x,y,_ = guide.project(xyz)
            draw.ellipse((x-4,y-4,x+4,y+4),fill=ORANGE,outline="white",width=1)
    for n,ident in enumerate(examples,1):
        x,y,_ = guide.project(np.asarray(base["parts"][ident]["box"]).mean(axis=0))
        xp,yp = x+22,max(25,y-26)
        draw.line((x,y,xp,yp),fill=INK,width=2)
        draw.ellipse((xp-13,yp-13,xp+13,yp+13),fill=INK,outline="white",width=2)
        label(draw,(xp,yp),str(n),16,"white",True,"mm")
    label(draw,(18,14),"整箱安装位置",21,bold=True)
    return view


def for_detail(ident, row, edges, parts):
    local = [edges[k] for k in row["evidence"] if k in edges]
    neighbors = []
    for edge in local:
        for other in (edge["a"],edge["b"],edge.get("host")):
            if isinstance(other,str) and other != ident and other in parts and "box" in parts[other] and other not in neighbors:
                neighbors.append(other)
    return neighbors[:3]


def detail(model, mesh, base, ident, neighbors, width, height, number):
    part = np.asarray(base["parts"][ident]["box"])
    center = part.mean(axis=0)
    scale = min(950,max(180,4*np.ptp(part,axis=0).max()))
    y = np.mean(np.asarray(base["frame_box"])[:,1])
    direction = [.45,.85 if center[1] >= y else -.85,.64]
    cx = Camera(width,height,center,direction,scale,up=(0,0,1))
    items = [ident,*neighbors]
    picture = render(cx,[mesh(i,.16 if i==ident else .65) for i in items],[ORANGE,BLUE,GRAY,GRAY][:len(items)])
    draw = ImageDraw.Draw(picture)
    draw.rectangle((0,0,width,35),fill="#edf3f7")
    label(draw,(12,7),f"{number}　{ident} 局部",18,bold=True)
    lines = []
    for n,other in enumerate(neighbors):
        try:
            witness = gap(model.local_shape(ident),model.local_shape(other))
        except (ValueError,RuntimeError,AssertionError):
            lines.append((other,None))
            continue
        p,q = cx.project([witness["point_1_mm"],witness["point_2_mm"]])[:,:2]
        if (0<=p[0]<width and 35<=p[1]<height-94):
            draw.line((*p,*q),fill="#673478",width=3)
            x,y = p
            draw.ellipse((x-5,y-5,x+5,y+5),outline="white",fill="#673478",width=1)
        lines.append((other,witness["distance_mm"]))
    draw.rectangle((0,height-90,width,height),fill="#f8fbfd")
    label(draw,(12,height-82),"橙色为待复核件，蓝灰色为相邻件",16)
    for n,(other,distance) in enumerate(lines):
        info = f"最短间隙 {distance:.3f} mm" if distance is not None else "量测失败，需核对"
        label(draw,(12,height-60+n*18),f"{other}　{info}",15,color=BLUE if n==0 else INK)
    if not lines:
        label(draw,(12,height-53),"暂无可用配对见证",16)
    return picture,lines


def render_review(analysis, rules, folder):
    package,_,source,raw = read_analysis(analysis)
    folder = Path(folder).resolve()
    folder.mkdir(parents=True)
    base = load(package/"base.json")
    groups = {rule:[i for i,r in raw["results"].items() if r["rule"]==rule and r["action"]=="暂留待判"] for rule in dict.fromkeys(rules)}
    groups = {rule:ids for rule,ids in groups.items() if ids}
    if not groups:
        save(folder/"截图索引.json",[])
        return
    model = Model(source)
    model.restore_basis(base)
    if model.paths.keys()!=raw["results"].keys():
        raise ValueError("Review CAD occurrences do not match the sealed prediction")
    by_edge = {e["edge_id"]:e for e in raw["interfaces"]}

    @lru_cache(maxsize=800)
    def mesh(ident, tolerance):
        return triangles(model.local_shape(ident),tolerance)

    frame = [mesh(i,4.) for i,r in raw["results"].items() if r["rule"] in FRAME_ROLES and "box" in base["parts"][i]]
    frame = [s for s in frame if len(s)]
    if not frame:
        raise ValueError("No frame geometry available to contextualize review images")
    receipts = []
    for rule,ids in groups.items():
        available = [i for i in ids if "box" in base["parts"][i]]
        if not available:
            receipts.append({"rule":rule,"count":len(ids),"image":None,"reason":"缺可绘制几何"})
            continue
        start = available[0]
        get_center = lambda i:np.asarray(base["parts"][i]["box"]).mean(axis=0)
        examples = [start]
        if len(available)>1:
            def host_role(ident):
                for edge_id in raw["results"][ident]["evidence"]:
                    host = by_edge.get(edge_id,{}).get("host")
                    if isinstance(host,str) and host in raw["results"]:
                        return raw["results"][host]["rule"]
                return None

            start_host = host_role(start)
            alternatives = [i for i in available[1:] if host_role(i)!=start_host] or available[1:]
            examples.append(max(alternatives,key=lambda i:np.linalg.norm(get_center(i)-get_center(start))))
        page = Image.new("RGB",(PAGE_WIDTH,900),"white")
        draw = ImageDraw.Draw(page)
        label(draw,(MARGIN,18),f"{raw['results'][start]['role']} · 待复核",31,bold=True)
        label(draw,(MARGIN,60),f"{len(ids)} 个候选，按源 STEP 展示安装位置。橙色点为待复核实例。",20)
        width = PAGE_WIDTH-2*MARGIN
        page.paste(overview(base,frame,available,examples,width,348),(MARGIN,101))
        column = (width-12)//len(examples)
        evidence = []
        for n,ident in enumerate(examples,1):
            neighbors = for_detail(ident,raw["results"][ident],by_edge,base["parts"])
            picture,lines = detail(model,mesh,base,ident,neighbors,column,375,n)
            page.paste(picture,(MARGIN+(n-1)*(column+12),462))
            evidence.append({"id":ident,"neighbors_and_gaps_mm":lines})
        label(draw,(MARGIN,854),"确认该组的真实连接与当前工况清理范围。",19)
        name = f"{rule}-待复核.png"
        page.save(folder/name)
        receipts.append({"rule":rule,"count":len(ids),"image":name,"examples":evidence})
    save(folder/"截图索引.json",receipts)
    fields = ["判据","匿名实例","角色","复核依据"]
    write_csv(folder/"截图覆盖实例.csv",fields,[dict(zip(fields,(rule,i,raw["results"][i]["role"],
              raw["results"][i]["falsifier"]))) for rule,ids in groups.items() for i in ids])
    print("Grouped review images",[(r["rule"],r["count"],r["image"]) for r in receipts],flush=True)
