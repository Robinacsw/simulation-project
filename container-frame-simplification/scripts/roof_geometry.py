"""Measure a cover-bracket and its local underside holder in box coordinates."""

from collections import defaultdict
from itertools import combinations
from math import pi

import numpy as np

from relations import Graph, intervals


TOUCH_MM = .10
LOCAL_SEARCH_MM = 3.0
AXIS_LINEAR_TOLERANCE_MM = .05


def _bore_sections(model, instance):
    """Transform analytic cylindrical spans into neutral box coordinates."""
    definition = model.by_path[model.paths[instance]]["definition"]
    cylinders = model.features[definition].get("cylinders", [])
    rotation, translation = model.pose(instance)
    axes = model.axes @ rotation
    origin = model.axes @ (translation - model.origin)
    for face in cylinders:
        if not face["inner"] or face["angle"] < 2.9:
            continue
        first, last = [axes @ np.asarray(face[k]) + origin for k in ("a", "b")]
        delta = abs(last - first)
        direction = int(np.argmax(delta))
        transverse = np.delete(delta, direction)
        if delta[direction] < 2.0 or np.max(transverse) > AXIS_LINEAR_TOLERANCE_MM:
            continue
        lo, hi = (first, last) if first[direction] < last[direction] else (last, first)
        yield (direction, round(face["r"], 4), *lo.round(3), *hi.round(3)), face["angle"]


def _round_bores(sections):
    pieces = defaultdict(list)
    for key, angle in sections:
        pieces[key].append(angle)
    for (direction, radius, ax, ay, az, bx, by, bz), angles in pieces.items():
        sweep = sum(angles)
        if abs(sweep - 2 * pi) <= .10:
            yield {"axis_start_mm":[ax, ay, az], "axis_end_mm":[bx, by, bz],
                   "direction":direction, "radius_mm":radius, "scope_degrees":round(sweep * 180 / pi, 4)}


def _vertical_slots(sections):
    halves = [key for key,angle in sections if key[0] == 2 and abs(angle - pi) <= .10]
    for first, second in combinations(halves, 2):
        _, radius, ax, ay, az, bx, by, bz = first
        _, other, cx, cy, cz, dx, dy, dz = second
        width = 2 * radius
        length = float(np.linalg.norm(np.asarray([cx-ax, cy-ay])))
        if (abs(other-radius) <= .03 and abs(cz-az) <= .03 and abs(dz-bz) <= .03
                and .65*width <= length <= 2.0*width):
            yield {"axis_start_mm":[ax, ay, az], "axis_end_mm":[bx, by, bz],
                   "second_axis_start_mm":[cx, cy, cz], "second_axis_end_mm":[dx, dy, dz],
                   "width_mm":width, "radius_mm":radius, "end_center_distance_mm":round(length,5)}


def _classify_shapes(graph):
    connectors, holders, covers, structural = {}, {}, [], []
    for i, ident in enumerate(graph.ids):
        if graph.center[i, 2] < graph.frame[0][2] + .87 * graph.span[2]:
            continue
        x, y, z = graph.size[i]
        row = graph.parts[ident]
        planar_min, planar_max = sorted((x, y))
        occupancy = row["volume"] / max(x * y * z, 1)
        if (.57 * graph.span[1] <= y <= .95 * graph.span[1]
                and .025 * graph.span[0] <= x <= .14 * graph.span[0]
                and 0 < z <= .10 * min(x, y) and occupancy <= .16
                and row["area"] > 0 and 2 * row["volume"] / row["area"] <= .32 * z):
            covers.append(i)
        if graph.types[ident] in {"TOP_RAIL", "CROSS"} and planar_max >= .48 * graph.span[1]:
            structural.append(i)
        if not (18 <= z <= 40 and 22 <= planar_min <= planar_max <= 80 and .16 <= occupancy <= .46):
            continue
        sections = list(_bore_sections(graph.model, ident))
        if .82 <= planar_min / planar_max <= 1.0:
            slots = [slot for slot in _vertical_slots(sections) if .10 <= slot["radius_mm"] / planar_min <= .20]
            if len(slots) == 1:
                connectors[i] = slots[0]
        if .49 <= planar_min / planar_max <= .82:
            bores = [bore for bore in _round_bores(sections) if bore["direction"] == 2
                     and .20 <= bore["radius_mm"] / planar_min <= .36]
            if len(bores) == 1:
                holders[i] = bores[0]
    return connectors, holders, np.asarray(covers, dtype=int), np.asarray(structural, dtype=int)


def _near_candidates(graph, index, options, tolerance):
    return options[np.max(intervals(graph.lo[options], graph.hi[options],
                                    graph.lo[index], graph.hi[index]), axis=1) <= tolerance]


def _cover_witnesses(graph, index, options):
    ident = graph.ids[index]
    evidence = []
    for other in _near_candidates(graph, index, options, TOUCH_MM):
        neighbor = graph.ids[other]
        if graph.parts[ident]["parent"] != graph.parts[neighbor]["parent"]:
            continue
        measured = graph.nearest(ident, neighbor)
        if measured is not None and measured["gap_mm"] <= TOUCH_MM and not measured["inner"]:
            evidence.append({"instance":neighbor, **measured})
    return evidence


def _host_witnesses(graph, index, options):
    ident = graph.ids[index]
    matches = []
    for other in _near_candidates(graph, index, options, LOCAL_SEARCH_MM):
        if graph.center[other, 2] >= graph.center[index, 2]:
            continue
        gap = graph.nearest(ident, graph.ids[other])
        if gap is not None and gap["gap_mm"] <= LOCAL_SEARCH_MM:
            matches.append((other, gap))
    return matches


def _slot_alignment(slot, bore):
    start, end = [np.asarray(slot[k][:2]) for k in ("axis_start_mm", "second_axis_start_mm")]
    center = np.asarray(bore["axis_start_mm"][:2])
    travel = end - start
    ratio = float((center-start) @ travel / (travel @ travel))
    offset = float(np.linalg.norm(center - (start + ratio * travel)))
    return ratio, offset


def _holder_witnesses(graph, index, options, slot, bores):
    ident = graph.ids[index]
    matches = []
    for other in _near_candidates(graph, index, options, LOCAL_SEARCH_MM):
        delta = graph.lo[index, 2] - graph.hi[other, 2]
        if not (-TOUCH_MM <= delta <= LOCAL_SEARCH_MM):
            continue
        overlap = np.minimum(graph.hi[index, :2], graph.hi[other, :2]) - np.maximum(graph.lo[index, :2], graph.lo[other, :2])
        if np.any(overlap < .60 * np.minimum(graph.size[index, :2], graph.size[other, :2])):
            continue
        holder = graph.ids[other]
        if graph.parts[ident]["parent"] == graph.parts[holder]["parent"]:
            continue
        ratio, offset = _slot_alignment(slot, bores[other])
        if not (-.06 <= ratio <= 1.06 and offset <= .15):
            continue
        gap = graph.nearest(ident, holder)
        if (gap is not None and gap["gap_mm"] <= LOCAL_SEARCH_MM and not gap["inner"]
                and abs(gap["gap_mm"] - max(0, delta)) <= TOUCH_MM):
            matches.append((other, gap, {"slot_axis_fraction":round(ratio,5), "offset_mm":round(offset,5)}))
    return matches


def measure_roof_geometry(model, base):
    """Use neutral assembly relationships and B-Rep surfaces for local roof fittings."""
    graph = Graph(model, base)
    connectors, holders, covers, structural = _classify_shapes(graph)
    holder_indices = np.asarray(list(holders), dtype=int)
    passed = []
    for index, bracket_hole in connectors.items():
        connector = graph.ids[index]
        cover_witnesses = _cover_witnesses(graph, index, covers)
        if not cover_witnesses:
            continue
        attachment = _holder_witnesses(graph, index, holder_indices, bracket_hole, holders)
        direct_hosts = _host_witnesses(graph, index, structural)
        if len(attachment) != 1 or direct_hosts:
            continue
        holder_index, connection, aligned = attachment[0]
        holder = graph.ids[holder_index]
        hosts = _host_witnesses(graph, holder_index, structural)
        if len(hosts) != 1:
            continue
        host_index, bearing = hosts[0]
        if bearing["gap_mm"] > TOUCH_MM or bearing["inner"]:
            continue
        passed.append({"connector":connector, "holder":holder, "host":graph.ids[host_index],
                       "cover_panels":[r["instance"] for r in cover_witnesses],
                       "bracket_hole":bracket_hole, "holder_hole":holders[holder_index], "slot_alignment":aligned,
                       "connector_to_holder":connection, "connector_to_covers":cover_witnesses,
                       "holder_to_host":bearing, "qualified":True,
                       "source":"B-Rep slot, round bore and nearest-point relationships"})
    by_holder = defaultdict(list)
    for relation in passed:
        by_holder[relation["holder"]].append(relation["connector"])
    result = []
    for evidence in passed:
        connector, holder = evidence["connector"], evidence["holder"]
        if len(by_holder[holder]) != 1:
            continue
        result.append(dict(evidence, a=connector, b=holder, kind="accessory_geometry", subtype="roof_connector"))
        result.append(dict(evidence, a=holder, b=connector, kind="accessory_geometry", subtype="roof_subassembly"))
    return result
