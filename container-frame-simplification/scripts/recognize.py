"""Frozen role and action rules for both anonymous modes. Standard library only."""

from collections import Counter, defaultdict
from math import dist, prod

RULE_VERSION = "frame-EU-1.0"

# Ratios from the previous anonymous geometry screen, and new pilot screens.
# These limit which candidates are examined. No rule here sets a weld tolerance.
PARAMETERS = {
    "corner_long_fraction": .10, "corner_width_fraction": .15,
    "end_layer_fraction": .16, "upright_height_fraction": .45,
    "member_aspect": 3.35, "body_distance_numeric_mm": .03,
    "port_search_mm": 9.0, "frame_neighbor_search_mm": 4.0,
    "same_corner_horizontal_mm": 105.0, "min_end_projection_mm2": 30.0,
    "linear_section_bbox_mm": 20.0, "long_rod_length_mm": 900.0,
    "guide_separation_fraction": .18, "short_joint_section_mm": 750.0,
}

ROLE_GROUPS = {
    "角件": "角件", "角柱": "框架立柱", "框架立柱候选": "框架立柱",
    "顶侧梁": "框架纵梁", "底侧梁": "框架纵梁", "底架纵梁": "框架纵梁",
    "门楣或顶端梁": "框架横梁", "端框底梁": "框架横梁", "底横梁": "框架横梁",
    "顶部横梁": "框架横梁", "分段框架梁": "框架横梁",
    "底板工作面板": "底板", "底封板": "底板",
    "斜撑": "节点加强", "框架跨接板": "节点加强",
    "门锁杆": "开闭附件", "门锁导向件": "开闭附件",
    "壳体板件": "围护板件", "局部圆柱连接件": "连接小件",
    "活动梁型连接组": "结构连接",
}


def snapshot(role, action, code, reason, falsifier, state):
    return {"role":role, "family":ROLE_GROUPS.get(role), "action":action,
            "rule":code, "basis":reason, "falsifier":falsifier, "state":state,
            "evidence":[]}


def coords(p, limits):
    low, high = p["box"]
    xyz = tuple((a+b)*.5 for a,b in zip(low,high))
    dims = tuple(b-a for a,b in zip(low,high))
    span = tuple(b-a for a,b in zip(*limits))
    pos = tuple((c-a)/r for c,a,r in zip(xyz,limits[0],span))
    return xyz,dims,span,pos


def casting_candidate(p, limits):
    if "box" not in p or p["solids"] != 1:
        return False
    _,_,_,(x,y,z) = coords(p,limits)
    a,b,c = p["shape_mm"]
    holes = p["opening_dirs"]
    multi_opening = any(sum(u*v for u,v in zip(m,n))<.30 for m in holes for n in holes)
    t = PARAMETERS
    return (min(x,1-x)<=t["corner_long_fraction"] and min(y,1-y)<=t["corner_width_fraction"]
            and min(z,1-z)<=t["end_layer_fraction"] and 65<=a<=b<=c<=320 and
            p["faces"]>=18 and p["volume"]>=430000 and multi_opening)


def upright_candidate(p, limits):
    if "box" not in p or p["solids"] != 1:
        return False
    _, box, span, (x,y,_) = coords(p,limits)
    a,b,c = p["shape_mm"]
    return (box[2]>=PARAMETERS["upright_height_fraction"]*span[2]
            and abs(p["direction"][2])>=.90 and c>=3.35*b and a>=16
            and (min(x,1-x)<.105 or min(y,1-y)<.19))


def geometry(p, limits):
    if "box" not in p:
        return snapshot("未判", "暂留待判", "READ", "实体特征读取失败", "补齐该实例及位姿", "读取失败")
    _, run, span, (x,y,z) = coords(p,limits)
    a,b,c = p["shape_mm"]
    edge_x, edge_y = min(x,1-x),min(y,1-y)
    slender = c >= PARAMETERS["member_aspect"] * max(1,b)
    t = 2*p["volume"]/max(1,p["area"])
    major = max(range(3),key=run.__getitem__)
    flat = a<=max(12,.04*b) and t<=max(10,.03*b)
    if p["solids"]!=1:
        return snapshot("多实体单元", "暂留待判", "MULTI", f"同一实例含 {p['solids']} 个实体；需要逐体分工或分割",
                        "分体核查可确定角色和处置", "需局部拆分")
    if not p["valid"]:
        return snapshot("未判", "暂留待判", "QUALITY", "B-Rep 有效性检查未通过", "复核原生几何和面积", "几何待核")
    opts = [
        (casting_candidate(p,limits), "角件","保留几何","CAST",
         "端侧上下交界的紧凑铸型，具有两个不平行的内圆弧开口", "开口与试验接口形状不符或同位变体"),
        (upright_candidate(p,limits), "框架立柱候选","保留几何","POST",
         "靠端或靠侧，长轴竖直、跨过箱高的主要区间", "两端不接固定框架且独立导向件随门扇安装"),
        (flat and b>=.15*span[1] and run[2]<=max(14, .05*b) and z<=.20,
         "底封板" if z<.06 else "底板工作面板","保留几何","DECK",
         "底部广展水平薄板，位于底端面或底架上方", "板面实为设备局部底盖，或工作面落在其他铺板"),
        (slender and a>=22 and abs(p["direction"][0])>=.94 and c>=520 and edge_y<.095 and z>=.81,
         "顶侧梁","保留几何","TOP_RAIL","顶部侧缘纵向型件，具有可测的连续截面", "端部落在设备或可移除围护件"),
        (slender and a>=22 and abs(p["direction"][0])>=.94 and c>=520 and edge_y<.095 and z<=.19,
         "底侧梁","保留几何","BOTTOM_RAIL","底部侧缘纵向型件，具有可测的连续截面", "端部落在可移除设备托架"),
        (slender and a>=28 and abs(p["direction"][0])>=.94 and c>=760 and z<=.17,
         "底架纵梁","保留几何","BOTTOM_LONG","底架高度内纵向型件", "独立导轨的两端无框架载荷接口"),
        (slender and a>=22 and abs(p["direction"][1])>=.94 and run[1]>=.24*span[1] and
         (z>=.76 or z<=.21 or edge_x<.065),
         ("门楣或顶端梁" if z>=.76 else "端框底梁") if edge_x<.065 else ("顶部横梁" if z>=.76 else "底横梁"),
         "保留几何","CROSS","横向截面构件位于顶或底架高度／端框", "端头位于门扇、设备或独立盖板总成"),
        (flat and b>=130 and c>=450 and (z>=.80 or edge_x<.105 or edge_y<.095),
         "壳体板件","清理几何","SKIN","外围板状或薄折面，单件截面支持围护候选",
         "板端跨接两侧承力梁，且构成门洞边界或框架节点"),
        (c<=105 and b<=54 and p["concentric_candidate"] and .19<z<.80,
         "局部圆柱连接件","清理几何","SMALL","中部小型孔轴件；形状支持辅助连接初判",
         "核到双侧承载母件、活动梁紧固或结构载荷接口"),
    ]
    for match,role,action,rule,basis,falsifier in opts:
        if match:
            return snapshot(role,action,rule,basis,falsifier,"已作初判")
    return snapshot("未判","暂留待判","UNCLEAR","单件形态或局部位置无法排除其他用途",
                    "补齐端部、宿主、载荷接口或运动关系", "证据不足")


def resolve(ident, p, results, edges, parts, limits):
    result = results[ident]
    if "box" not in p or p["solids"]!=1 or not p["valid"]:
        return result
    ends = [e for e in edges if e["kind"]=="port" and e["a"]==ident and e["qualified"]]
    lows = [e for e in ends if e["side"]==0]
    highs = [e for e in ends if e["side"]==1]
    compatible = [(lo,hi) for lo in lows for hi in highs if
                  dist(lo["target_center"][:2],hi["target_center"][:2])<=PARAMETERS["same_corner_horizontal_mm"]]
    guide = [e for e in edges if e["kind"]=="axis_guide" and e["a"]==ident and e["panel"]]
    if compatible:
        lo,hi = max(compatible, key=lambda pair: (sum(e["area_mm2"] for e in pair), -sum(e["gap_mm"] for e in pair)))
        pick = [lo,hi]
        possible_pairs = sorted({(e["b"],f["b"]) for e,f in compatible})
        out = snapshot("角柱","保留几何","POST_TWO",
                       "上下端有效截面分别投影到同角位的两角件，几何定位已核",
                       "一端失去端面投影，或归入不同纵向角位", "两端配对已量测")
        out["station_duplicates"] = [list(pair) for pair in possible_pairs] if len(possible_pairs)>1 else []
        if guide:
            out["action"]="暂留待判"
            out["state"]="连接线索冲突"
            pick += guide
        out["evidence"] = [e["edge_id"] for e in pick]
        return out
    if guide:
        by_panel=defaultdict(list)
        for g in guide:
            by_panel[g["panel"]].append(g)
        pairs=[(lo,hi) for links in by_panel.values() for lo in links for hi in links if
               hi["axis_height"]-lo["axis_height"]>=PARAMETERS["guide_separation_fraction"]*p["shape_mm"][2]]
        if pairs:
            low,high = max(pairs,key=lambda es:es[1]["axis_height"]-es[0]["axis_height"])
            out = snapshot("门锁杆","清理几何","LOCK_GUIDE","两个分开的导向界面均与杆轴同心搭接，且同一门板型宿主可定位",
                           "截面直接通到角件，或导向件和门板实际分属不同机构","几何关系已量测")
            out["evidence"]=[low["edge_id"],high["edge_id"]]
            if ends:
                out["action"]="暂留待判"
                out["state"]="端面与导向线索冲突"
                out["evidence"] += [e["edge_id"] for e in ends]
            return out
    if len(ends) and result["rule"]=="POST":
        out = snapshot("框架立柱候选","暂留待判","POST_ONE",
                       "角位立件已有一端配对；另一端或同角位条件需核实",
                       "补齐另一端实体配对，或查得门扇宿主", "关键端部证据不足")
        out["evidence"]=[e["edge_id"] for e in ends]
        return out
    joints = [e for e in edges if e["kind"]=="beam_clamp" and e["a"]==ident]
    if joints:
        es=[(b,f) for b in joints for f in joints if b["b"]!=f["b"] and b["long_end"] and
            b["shaft_cylinder"]==f["shaft_cylinder"]]
        if es:
            beam,side=max(es,key=lambda pair:(min(e["axial_overlap"] for e in pair),-sum(e["center_offset"] for e in pair)))
            out=snapshot("活动梁型连接组","保留几何","CLAMP_TWO",
                         "杆轴分别进入梁端和另一母件的内圆柱面；保留候选连接及两侧母件",
                         "两侧实际属于同一母件、螺杆未夹持或梁与外侧件可永久连接","夹持关系候选，实际活动自由度待核")
            out["evidence"]=[beam["edge_id"],side["edge_id"]]
            out["unit_hosts"]=[beam["b"],side["b"]]
            return out
        if result["action"]!="保留几何":
            out=snapshot("局部圆柱连接件","暂留待判","CLAMP_ONE",
                         "有梁端或接板的孔轴配对，须核到另一母件才能判结构连接","补齐另一侧夹持件","连接证据不足")
            out["evidence"]=[e["edge_id"] for e in joints]
            return out
    if any(e["kind"]=="axis_guide" and e["b"]==ident and e["panel"] and
           results[e["a"]]["rule"]=="LOCK_GUIDE" and results[e["a"]]["action"]=="清理几何"
           for e in edges):
        used=[e["edge_id"] for e in edges if e["kind"]=="axis_guide" and e["b"]==ident
              and e["panel"] and results[e["a"]]["rule"]=="LOCK_GUIDE"]
        out=snapshot("门锁导向件","清理几何","GUIDE","和已识别门锁杆同轴，导向位置可定位",
                     "支座也连接框架焊接节点或活动梁","几何关系已量测")
        out["evidence"]=used
        return out
    host_edges = [e for e in edges if e["kind"]=="frame_end" and e["a"]==ident and
                  e["gap_mm"]<=PARAMETERS["body_distance_numeric_mm"] and not e["inner"]]
    low_end = [e for e in host_edges if e["side"]==0]
    high_end = [e for e in host_edges if e["side"]==1]
    if not low_end or not high_end:
        return result
    selected = next(((a,b) for a in low_end for b in high_end if a["b"]!=b["b"]),None)
    if selected is None:
        return result
    _,box,span,(x,y,z) = coords(p,limits)
    a,b,c = p["shape_mm"]
    if result["rule"]=="SKIN" and z>=.76 and box[1]>=.43*span[1] and box[0]<.23*span[1]:
        out=snapshot("框架跨接板","保留几何","BRIDGE",
                     "短进深的顶横向板两边分别接近不同的顶侧纵梁；保留以核门楣及节点作用",
                     "两端母件是可拆盖板，或板边未焊接且只作围护","几何两端已量测，工艺待核")
    elif a>=19 and c>=3.35*max(b,1) and result["rule"]=="UNCLEAR":
        out=snapshot("分段框架梁","保留几何","FRAME_SEG",
                     "轴向两个端部各定位到不同框架母件，保留端部传力截面",
                     "端部为非承载卡扣，或测得超出节点的接近","几何两端已量测，工艺待核")
    else:
        return result
    out["evidence"]=[e["edge_id"] for e in selected]
    return out


def predict(packet, prompt_sha):
    mode = packet["mode"]
    data = packet["base"]
    assert mode in ("G","G+C") and packet["rules"] == RULE_VERSION
    assert packet["authorization"]=="只读本任务包"
    assert set(packet)=={"mode","base","rules","authorization","cases","scope","relations","run"}
    assert set(data)=={"model","configuration","unit","basis","frame_box","parts","tree"}
    assert data["unit"]=="mm" and data["configuration"]=="C01"
    assert all(set(p)<={"parent","availability","box","shape_mm","direction","solids","faces",
                       "volume","area","valid","shaft","opening_dirs","concentric_candidate"}
               for p in data["parts"].values())
    assert all(ident.startswith("I") and ident[1:].isdigit() for ident in data["parts"])
    assert all(ident.startswith("T") and ident[1:].isdigit() for ident in data["tree"])
    edges = packet["relations"]
    assert (mode=="G+C" or edges is None) and (mode=="G" or isinstance(edges,list))
    results={ident:geometry(p,data["frame_box"]) for ident,p in data["parts"].items()}
    if mode=="G":
        for row in results.values():
            row["relation_state"]="本模式有意屏蔽"
    else:
        adjacency=defaultdict(list)
        for edge in edges:
            adjacency[edge["a"]].append(edge)
            adjacency[edge["b"]].append(edge)
        for ident,p in data["parts"].items():
            if ident in adjacency:
                results[ident] = resolve(ident,p,results,adjacency[ident],data["parts"],data["frame_box"])
        # Resolve guides after all rod roles are available, independent of ID order.
        for ident,p in data["parts"].items():
            if ident in adjacency and results[ident]["rule"] not in {"LOCK_GUIDE","POST_TWO","POST_ONE","BRIDGE","FRAME_SEG","CLAMP_TWO"}:
                results[ident] = resolve(ident,p,results,adjacency[ident],data["parts"],data["frame_box"])
        for ident,row in results.items():
            used = row["evidence"]
            row["relation_state"]="几何有证据，施工及仿真连接待验证" if used else (
                "有测量记录，未形成有效配对" if ident in adjacency else "本轮未取得关系证据")
    groups=[]
    for ident,result in results.items():
        groups.append({"eu": f"U{ident[1:]}", "role": result["role"], "members":[ident],
                       "interfaces": result["evidence"], "action":result["action"]})
    if mode=="G+C":
        for ident,result in results.items():
            if result["rule"]=="LOCK_GUIDE":
                selected = [e for e in adjacency[ident] if e["edge_id"] in result["evidence"]]
                groups.append({"eu":f"D{ident[1:]}","role":"门锁杆导向机构",
                               "members":sorted({ident}|{e["b"] for e in selected if e["kind"]=="axis_guide"}),
                               "interfaces": result["evidence"],"action":result["action"]})
            if result["rule"]=="CLAMP_TWO":
                groups.append({"eu":f"J{ident[1:]}","role":"梁端孔轴夹持候选",
                               "members":[ident,*result["unit_hosts"]],"interfaces":result["evidence"],
                               "action":result["action"]})
    return {"run":packet["run"],"model":data["model"],"mode":mode,"prompt_sha256":prompt_sha,
            "rules":RULE_VERSION, "results":results, "units":groups,
            "interfaces":edges if mode=="G+C" else [],
            "actions":Counter(r["action"] for r in results.values())}
