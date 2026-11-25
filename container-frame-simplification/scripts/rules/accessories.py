"""Rules appended to the frozen G+C recognizer; no engineering-answer access."""

accessory_predecessor = predict
ACCESSORY_REVISION = 'frame-EU-1.3'

ACCESSORY_ROLE = {
    'end_corrugation': ('端部波板', 'END_CORRUGATION', '箱体端平面的大幅薄折面，周期性竖向折弯与同层分幅边界'),
    'alarm': ('报警器类外设', 'ALARM_PACKAGE', '局部多实体外壳与内部小件，未形成框架连续边或母件夹持链'),
    'mesh': ('钢网', 'WIRE_MESH', '重复细长网丝构成两组近正交行列，分布于围护平面'),
    'clip': ('绑线扣候选', 'CABLE_CLIP', '带通孔的低体量双脚薄带，架空段可走线，作用域限于局部'),
    'ring': ('系固绳环候选', 'LASH_RING', '连续弯环有两个共面脚，局部固定，未形成框架跨接'),
    'side_hinge': ('独立侧门铰链座', 'SIDE_HINGE', '沿长侧门线的两耳竖向孔与中间转轴和转耳同轴，且与端门位置区分'),
    'bottom_pad': ('底架垫板', 'BOTTOM_PAD', '单根底侧梁下的平面贴合板，按底角件离散支承边界判读'),
    'fixture_angle': ('侧固定角钢', 'FIXTURE_ANGLE', '底侧梁外侧 L 截面与局部两片加劲组成单母件组'),
    'fixture_gusset': ('固定角钢加强片', 'FIXTURE_GUSSET', '三角加劲片的贴合面闭合于同一外侧角钢组'),
}
FIXTURES = {'bottom_pad', 'fixture_angle', 'fixture_gusset'}
def accessory_disposition(evidence, old, packet):
    if old['rule'] in {'CLAMP_REVIEW', 'PANEL_NEEDS_ATTACHMENT'}:
        return '暂留待判', '已有框架、跨接或夹持证据与附件签名冲突'
    if len(evidence) != 1:
        return '暂留待判', '出现多个附件用途候选'
    edge = evidence[0]
    if edge.get('qualified') is False:
        return '暂留待判', '形状吻合，宿主或邻件证据未闭合'
    if edge['subtype'] in {'mesh','ring','clip'}:
        if not any(w['gap_mm'] <= .10 for w in edge['neighbor_witnesses']):
            return '暂留待判', '已定位附件形状，尚无数值贴合的宿主或网丝见证'
        if edge['subtype'] in {'ring','clip'} and sum(w['gap_mm'] <= .10 for w in edge['neighbor_witnesses']) != 1:
            return '暂留待判', '多个接近的母件竞争，检查局部夹持或跨接'
    if edge['subtype'] in FIXTURES:
        if edge.get('weld_interface_review'):
            return '暂留待判', '发现与另一保留构件相对面及小间隙，核对是否为预留焊接界面'
        if packet['scope'].get('fixture_boundary') != 'discrete_castings_and_floor' or not edge.get('host'):
            return '暂留待判', '需先核对支承、限位、吊点和夹具作用面'
    if edge['subtype'] == 'alarm':
        return '清理几何', '端侧设备包证据通过，局部壳体抽测用于邻接对账'
    return '清理几何', '形状、箱体局部位置及附件关系通过'


def predict(packet, prompt_sha):
    assert packet['rules'] == ACCESSORY_REVISION
    evidence = defaultdict(list)
    previous_edges = []
    for edge in packet['relations']:
        if edge['kind'] == 'accessory_geometry':
            assert edge['subtype'] in ACCESSORY_ROLE
            evidence[edge['a']].append(edge)
        else:
            previous_edges.append(edge)
    result = accessory_predecessor(dict(packet, rules=FINAL_REVISION, relations=previous_edges), prompt_sha)
    values = result['results']
    for ident, relations in evidence.items():
        role, code, basis = ACCESSORY_ROLE[relations[0]['subtype']]
        old = values[ident]
        action, state = accessory_disposition(relations, old, packet)
        # An existing keep remains in force while a conflicting new signature is reviewed.
        if old['action'] == '保留几何':
            old['accessory_conflict'] = [edge['edge_id'] for edge in relations]
            continue
        values[ident] = snapshot(role, action, code, basis,
                                 '跨两主构件、落入承载面、存在结构夹持或用途竞争时必须复核', state)
        values[ident].update(family='局部附件或围护', evidence=[edge['edge_id'] for edge in relations],
                             relation_state=state)
    for eu in result['units']:
        if len(eu['members']) == 1:
            row = values[eu['members'][0]]
            eu.update(role=row['role'], action=row['action'], interfaces=row['evidence'])
    result['interfaces'] = packet['relations']
    result['rules'] = ACCESSORY_REVISION
    result['actions'] = Counter(row['action'] for row in values.values())
    return result
