"""Appended to the frozen 1.0 rules in the isolated G+C process."""

prior_predict = predict
REVISION = "frame-EU-1.1"


def measured_stack(edge, parts):
    hosts=edge["hosts"]
    links=edge.get("links",[])
    fasteners=edge["fasteners"]
    if not edge["qualified"] or len(hosts)!=2 or len(fasteners)<3 or edge.get("other_host_pairs",0)>0:
        return False
    plates=[h["instance"] for h in hosts]
    backs=[h["backing"]["parts"][1] for h in hosts]
    if len(set(plates+backs))!=4 or set(plates+backs)&set(fasteners):
        return False
    if set(fasteners+plates+backs)-parts.keys() or len(set(fasteners))!=len(fasteners):
        return False
    if any(h["axial_overlap_mm"]<1.0 or len(h["mouth_faces"])!=2 or
           h.get("shaft_probe_common_volume_mm3",1)>.05 or
           h["backing"]["parts"][0]!=h["instance"] or
           (h["backing"].get("projected_area_mm2",0)<12 and h["backing"].get("common_volume_mm3",0)<=1)
           for h in hosts):
        return False
    if len(links)!=len(fasteners)+1 or links[0]["parts"][0]!=edge["a"]:
        return False
    if any(pair["projected_area_mm2"]<=0 for pair in links):
        return False
    if any(a["parts"][1]!=b["parts"][0] for a,b in zip(links,links[1:])):
        return False
    layers=[links[0]["parts"][0],*[edge["parts"][1] for edge in links]]
    return len(set(layers))==len(layers) and set(layers)==set(fasteners+plates)


def predict(packet, prompt_sha):
    assert packet["mode"]=="G+C" and packet["rules"]==REVISION
    result=prior_predict(dict(packet,rules=RULE_VERSION),prompt_sha)
    values=result["results"]
    data=packet["base"]["parts"]
    rels=packet["relations"]
    active=defaultdict(list)
    incomplete=defaultdict(list)
    host_of=defaultdict(list)
    thin=defaultdict(list)
    for e in rels:
        if e["kind"]=="shaft_stack":
            if measured_stack(e,data):
                for i in e["fasteners"]:
                    active[i].append(e)
                for side in e["hosts"]:
                    for i in side["backing"]["parts"]:
                        host_of[i].append(e["edge_id"])
            else:
                for i in e["fasteners"]:
                    incomplete[i].append(e["edge_id"])
        if e["kind"]=="panel_interface":
            thin[e["a"]].append(e["edge_id"])
    for i, evidence in active.items():
        ambiguous=len(evidence)>1
        plane_contact=all(e["same_plane"] for e in evidence)
        values[i]=snapshot("活动梁连接件",
                           "暂留待判" if ambiguous else "保留几何",
                           "CLAMP_AMBIGUOUS" if ambiguous else "CLAMP_STACK",
                           "同一螺杆穿过两块端板的孔，两侧端板均有独立长梁母件，已核螺栓头至螺母的有序面投影",
                           "缺一侧母件、孔不贯通、承压面失去交叠或多个套件竞争",
                           ("同一实例落入多套几何夹持" if ambiguous else
                            "作用面均在CAD数值贴合范围，预紧待核" if plane_contact else
                            "夹持层序已定位，CAD存在轴向间距或重叠"))
        values[i]["family"]="结构连接"
        values[i]["evidence"]=[e["edge_id"] for e in evidence]
        values[i]["relation_state"]="穿孔与承压面几何配对，工程夹紧状态待核"
    for i, evidence in host_of.items():
        if values[i]["action"]!="保留几何":
            values[i]=snapshot("夹持连接母件","保留几何","CLAMP_HOST",
                               "本件参与可定位的两侧端板与梁件接合",
                               "长梁母件重合或夹持链不成立","母件真实接合面已定位")
            values[i]["evidence"]=sorted(set(evidence))
            values[i]["relation_state"]="母件几何定位"
    for i, evidence in incomplete.items():
        if i in values and values[i]["action"]=="清理几何":
            values[i]=snapshot("夹持连接候选","暂留待判","CLAMP_REVIEW",
                               "发现穿孔及部分母件，完整两侧承压链待核",
                               "补齐两侧母件或证实为附件固定","部分几何证据")
            values[i]["evidence"]=sorted(set(evidence))
            values[i]["relation_state"]="母件或承压接口未配齐"
    for i, old in list(values.items()):
        if old["rule"]!="BRIDGE":
            continue
        # A thin panel touching beam ends still needs a documented attachment
        # or a physically identified mechanical connection to promote it.
        values[i]=snapshot("薄板跨接","暂留待判","PANEL_NEEDS_ATTACHMENT",
                           "两端最短距离或平面贴合，尚缺可追溯焊接或机械连接证据",
                           "核到真实承力焊缝、机械夹持或工况指定接口","构件作用待工程复核")
        values[i]["evidence"]=[*old["evidence"],*thin[i]]
        values[i]["relation_state"]="几何邻接，连接工艺待核"
    for eu in result["units"]:
        if len(eu["members"])==1:
            row=values[eu["members"][0]]
            eu.update(role=row["role"],interfaces=row["evidence"],action=row["action"])
    result["units"] += [
        {"eu":"V"+e["a"][1:],"role":"两侧长梁夹持连接",
         "members":sorted(set(e["fasteners"]+[i for h in e["hosts"] for i in h["backing"]["parts"]])),
         "interfaces":[e["edge_id"]],"action":"保留几何"}
        for e in rels if e["kind"]=="shaft_stack" and measured_stack(e,data)
        and all(values[i]["rule"]=="CLAMP_STACK" for i in e["fasteners"])]
    result["rules"]=REVISION
    result["actions"]=Counter(row["action"] for row in values.values())
    return result
