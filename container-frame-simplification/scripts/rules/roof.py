"""Anonymous project-role rules for roof mounts and approved bottom pads."""

roof_predecessor = predict
ROOF_REVISION = 'frame-EU-1.4'
ROOF_SUBTYPES = {'roof_connector', 'roof_subassembly'}


def confirmed_bottom_pad(edge, scope):
    if scope.get('fixture_boundary') != 'discrete_castings_and_floor' or not scope.get('bottom_pad_underskin_approved'):
        return False
    touching = edge['patch']
    if not edge.get('host') or edge['nearest']['inner'] or abs(touching['gap_mm']) > .03:
        return False
    # The authorized exception is a small lateral projection to the bottom skin.
    facing = [n for n in edge['additional_frame_neighbors'] if n['patch']]
    return (bool(facing) and all(n['role_shape']=='DECK' and abs(n['patch']['gap_mm']-3.0) <= .05
                                 and not n['nearest']['inner'] for n in facing)
            and sum(n['patch']['projected_area_mm2'] for n in facing) <= .03*touching['projected_area_mm2'])


def roof_evidence_ok(edge, reciprocal):
    if not (edge.get('qualified') and reciprocal.get('qualified')):
        return False
    if {edge['subtype'],reciprocal['subtype']} != ROOF_SUBTYPES:
        return False
    if edge['a'] != reciprocal['b'] or edge['b'] != reciprocal['a']:
        return False
    if len({edge['connector'],edge['holder'],edge['host']}) != 3:
        return False
    if edge['host'] != reciprocal['host']:
        return False
    mount = edge['connector_to_holder']
    host = edge['holder_to_host']
    return (.05 < mount['gap_mm'] <= 2.0 and not mount['inner'] and host['gap_mm'] <= .03 and not host['inner']
            and any(w['gap_mm'] <= .03 and not w['inner'] for w in edge['connector_to_covers']))


def predict(packet, prompt_sha):
    assert packet['rules']==ROOF_REVISION
    roof = defaultdict(list)
    inherited = []
    for edge in packet['relations']:
        if edge['kind']=='accessory_geometry' and edge['subtype'] in ROOF_SUBTYPES:
            roof[edge['a']].append(edge)
        else:
            inherited.append(edge)
    result = roof_predecessor(dict(packet, rules=ACCESSORY_REVISION, relations=inherited),prompt_sha)
    values = result['results']
    scope = packet['scope']
    for edge in inherited:
        if edge.get('subtype')!='bottom_pad' or not confirmed_bottom_pad(edge,scope):
            continue
        row = values[edge['a']]
        if row['rule']=='BOTTOM_PAD' and row['action']=='暂留待判':
            why = '单底侧梁贴合；仅邻近下层底蒙板且面投影较小；当前离散角件工况已确认清理'
            row.update(action='清理几何',state='已作初判',basis=why,relation_state=why,
                       falsifier='改用连续托底、出现其他母件接触、或超出已确认的底蒙板界面时复核')

    for ident, edges in roof.items():
        edge = edges[0]
        candidates = roof.get(edge['b'],[])
        confirmed = len(edges)==len(candidates)==1 and roof_evidence_ok(edge,candidates[0])
        row = values[ident]
        if row['action']=='保留几何' or row['rule'] in {'CLAMP_REVIEW','PANEL_NEEDS_ATTACHMENT'}:
            row['roof_accessory_conflict'] = [e['edge_id'] for e in edges]
            continue
        is_connector = edge['subtype']=='roof_connector'
        role = '顶盖连接板' if is_connector else '顶盖连接板下方固定座'
        code = 'ROOF_CONNECTOR' if is_connector else 'ROOF_HOLDER'
        basis = ('上方折弯带孔件贴薄顶盖，下方短固定座间隙已测；固定座只有一根框架母件'
                 if confirmed else '已取得屋顶安装件形状，配对面或单母件证据待核')
        adopted = is_connector and confirmed and scope.get('roof_connector_approved')
        action = '清理几何' if adopted else '暂留待判'
        state = ('按项目范围清理顶盖连接板' if adopted else
                 '配对固定座已定位，等待工程确认清理范围' if confirmed else '邻件或连接证据不足')
        values[ident] = snapshot(role,action,code,basis,'主梁或角件承担试验接口，出现第二个母件或配对缺失时复核',state)
        values[ident].update(family='可拆顶盖附件',evidence=[e['edge_id'] for e in edges],relation_state=state)

    for eu in result['units']:
        if len(eu['members'])==1:
            row = values[eu['members'][0]]
            eu.update(role=row['role'],action=row['action'],interfaces=row['evidence'])
    result['rules'] = ROOF_REVISION
    result['interfaces'] = packet['relations']
    result['actions'] = Counter(row['action'] for row in values.values())
    return result
