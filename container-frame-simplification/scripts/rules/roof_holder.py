"""Adopt verified local roof holders within the project cleanup scope."""

holder_predecessor = predict
HOLDER_REVISION = 'frame-EU-1.5'


def approved_local_holder(edge, reverse, results):
    holder, connector = edge['a'], edge['b']
    return (edge['subtype'] == 'roof_subassembly' and reverse['subtype'] == 'roof_connector'
            and edge['holder'] == reverse['holder'] == holder
            and edge['connector'] == reverse['connector'] == connector
            and roof_evidence_ok(edge, reverse) and roof_evidence_ok(reverse, edge)
            and results[holder]['rule'] == 'ROOF_HOLDER' and results[holder]['action'] == '暂留待判'
            and results[connector]['rule'] == 'ROOF_CONNECTOR' and results[connector]['action'] == '清理几何')


def predict(packet, prompt_sha):
    assert packet['rules'] == HOLDER_REVISION
    result = holder_predecessor(dict(packet, rules=ROOF_REVISION), prompt_sha)
    scope = packet['scope']
    result['rules'] = HOLDER_REVISION
    if not (scope.get('local_roof_holder_approved') is True
            and scope.get('fixture_boundary') == 'discrete_castings_and_floor'):
        return result

    roles = defaultdict(list)
    for edge in packet['relations']:
        if edge['kind'] == 'accessory_geometry' and edge['subtype'] in ROOF_SUBTYPES:
            roles[edge['a']].append(edge)

    values = result['results']
    adopted = set()
    for ident, edges in roles.items():
        if len(edges) != 1:
            continue
        edge = edges[0]
        reverse = roles.get(edge['b'], [])
        if len(reverse) != 1 or not approved_local_holder(edge, reverse[0], values):
            continue
        why = '与已清理的顶盖连接板有唯一槽轴配对，且仅贴合一根保留顶框梁；按本项目离散角件工况清理固定座'
        values[ident].update(action='清理几何', state='已作初判', basis=why, relation_state=why)
        adopted.add(ident)

    for unit in result['units']:
        if len(unit['members']) == 1 and unit['members'][0] in adopted:
            unit['action'] = values[unit['members'][0]]['action']
    result['actions'] = Counter(row['action'] for row in values.values())
    return result
