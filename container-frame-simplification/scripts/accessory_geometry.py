"""Measure anonymous accessory candidates against frozen B-Rep occurrences."""

from collections import Counter
from itertools import combinations

import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.GeomAbs import GeomAbs_BSplineSurface, GeomAbs_Plane
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
from OCP.TopoDS import TopoDS

from intake import indexed, properties
from measure_geometry import bounds, gap


def _rounded(values):
    return np.round(values, 5).tolist()


def _mesh(shape, box):
    thin_axis = int(np.argmin(np.ptp(box, axis=0)))
    rods = []
    solids = indexed(shape, TopAbs_SOLID)
    for solid in solids:
        rod_box = bounds(solid)
        sizes = np.ptp(rod_box, axis=0)
        axis = int(sizes.argmax())
        minor = np.delete(sizes, axis)
        if axis == thin_axis or sizes[axis] < 15 * max(minor):
            continue
        if min(minor) <= 0 or max(minor) > 1.10 * min(minor):
            continue
        length = sizes[axis]
        if length < .80 * (box[1, axis] - box[0, axis]):
            continue
        diameter = float(minor.mean())
        fill = abs(properties(solid, True)[0]) / (np.pi * (diameter / 2)**2 * length)
        if not .90 <= fill <= 1.10:
            continue
        transverse = ({0, 1, 2} - {thin_axis, axis}).pop()
        rods.append((axis, diameter, rod_box.mean(axis=0)[transverse]))
    directions = Counter(axis for axis, _, _ in rods)
    if len(rods) < .90 * len(solids) or len(directions) != 2 or min(directions.values()) < 4:
        return None
    diameters = [diameter for _, diameter, _ in rods]
    if max(diameters) > 1.20 * min(diameters):
        return None
    pitches = {}
    for axis in directions:
        centers = sorted(track for group, _, track in rods if group == axis)
        steps = np.diff(centers)
        pitch = float(np.median(steps))
        if pitch <= max(diameters) or np.count_nonzero(np.abs(steps-pitch) < .12*pitch) < .8*len(steps):
            return None
        pitches[str(axis)] = round(pitch, 5)
    return {
        "rod_count": len(rods),
        "rod_groups_by_local_axis": dict(directions),
        "median_rod_diameter_mm": round(float(np.median(diameters)), 5),
        "rod_pitch_by_local_axis_mm": pitches,
    }


def _faces(shape):
    for source in indexed(shape, TopAbs_FACE):
        face = TopoDS.Face_s(source)
        surface = BRepAdaptor_Surface(face, True)
        area, center = properties(face)
        normal_axis = None
        if surface.GetType() == GeomAbs_Plane:
            normal = np.abs(surface.Plane().Axis().Direction().Coord())
            if max(normal) > .99999:
                normal_axis = int(normal.argmax())
        yield face, surface.GetType(), area, np.array(center), normal_axis


def _paired_feet(faces, box, span_axis, normal_axis=None):
    dimensions = np.ptp(box, axis=0)
    center = box.mean(axis=0)
    feet = []
    for face, _, area, mid, axis in faces:
        if axis is None or axis == span_axis or area <= 10:
            continue
        if normal_axis is not None and axis != normal_axis:
            continue
        if min(abs(box[:, axis]-mid[axis])) > .10:
            continue
        if abs(mid[span_axis]-center[span_axis]) < .25*dimensions[span_axis]:
            continue
        if np.ptp(bounds(face), axis=0)[span_axis] > .30*dimensions[span_axis]:
            continue
        feet.append((axis, area, mid))
    pairs = []
    for (axis, area_a, a), (other_axis, area_b, b) in combinations(feet, 2):
        if axis != other_axis or abs(a[axis]-b[axis]) > .10:
            continue
        if abs(a[span_axis]-b[span_axis]) < .60*dimensions[span_axis]:
            continue
        if min(area_a, area_b) < .72 * max(area_a, area_b):
            continue
        pairs.append((area_a+area_b, axis, [area_a, area_b], np.array([a, b])))
    return max(pairs, key=lambda pair: pair[0], default=None)


def _ring(shape, box, data):
    sizes = np.ptp(box, axis=0)
    long_axis = int(sizes.argmax())
    faces = list(_faces(shape))
    dominant_curves = [
        area for face, kind, area, _, _ in faces
        if kind == GeomAbs_BSplineSurface
        and np.ptp(bounds(face), axis=0)[long_axis] > .70*sizes[long_axis]
    ]
    curved_ratio = sum(dominant_curves) / data["area"]
    if len(dominant_curves) < 2 or not .55 <= curved_ratio < .92:
        return None
    feet = _paired_feet(faces, box, long_axis)
    if feet is None:
        return None
    _, axis, areas, centers = feet
    return {
        "curved_span_surface_count": len(dominant_curves),
        "curved_span_surface_fraction": round(curved_ratio, 5),
        "mounting_face_axis_prototype": axis,
        "mounting_face_areas_mm2": _rounded(areas),
        "mounting_face_centers_prototype_mm": _rounded(centers),
    }


def _clip(shape, box, data):
    sizes = np.ptp(box, axis=0)
    center = box.mean(axis=0)
    long_axis = int(sizes.argmax())
    by_axis = {}
    for cylinder in data["cylinders"]:
        if not cylinder["inner"]:
            continue
        start, end = np.array(cylinder["a"]), np.array(cylinder["b"])
        delta = abs(end-start)
        axis = int(delta.argmax())
        mid = (start+end)/2
        sheet = delta[axis]
        if axis == long_axis or np.linalg.norm(np.delete(delta, axis)) > .001:
            continue
        if not .5 <= sheet <= .30*sizes[axis]:
            continue
        span_axes = [index for index in range(3) if index != axis]
        if np.any(abs(mid[span_axes]-center[span_axes]) > .13*sizes[span_axes]):
            continue
        diameter = 2*cylinder["r"]
        if not .20 <= diameter / min(sizes[span_axes]) <= .85:
            continue
        key = (axis, tuple(np.round(mid, 3)), round(diameter, 3), round(float(sheet), 3))
        by_axis[key] = by_axis.get(key, 0.0) + cylinder["angle"]
    holes = [key for key, angle in by_axis.items() if 6.1 <= angle <= 6.45]
    if len(holes) != 1:
        return None
    axis, center, diameter, sheet = holes[0]
    feet = _paired_feet(_faces(shape), box, long_axis, normal_axis=axis)
    if feet is None:
        return None
    _, _, areas, centers = feet
    clearance = abs(center[axis]-centers[0, axis])
    if clearance < 3*sheet or clearance < .55*sizes[axis]:
        return None
    return {
        "sheet_gauge_mm": sheet,
        "central_hole_diameter_mm": diameter,
        "central_hole_to_foot_plane_mm": round(float(clearance), 5),
        "mounting_face_axis_prototype": axis,
        "mounting_face_areas_mm2": _rounded(areas),
        "mounting_face_centers_prototype_mm": _rounded(centers),
    }


def _alarm_package(shape, box):
    size = np.ptp(box, axis=0)
    package_parts = []
    for solid in indexed(shape, TopAbs_SOLID):
        limits = bounds(solid)
        dimensions = np.ptp(limits, axis=0)
        volume = abs(properties(solid, True)[0])
        package_parts.append((dimensions, limits, volume))
    small = sum(np.all(dimensions <= .22*size) for dimensions, _, _ in package_parts)
    if small < .40*len(package_parts):
        return None
    for axis in range(3):
        plan = [dim for dim in range(3) if dim != axis]
        enclosing = [
            (dimensions, limits, volume)
            for dimensions, limits, volume in package_parts
            if np.all(dimensions[plan] >= .75*size[plan])
            and .25*size[axis] <= dimensions[axis] <= .80*size[axis]
            and volume < .35*np.prod(dimensions)
        ]
        for first, second in combinations(enclosing, 2):
            one, two = first[1], second[1]
            overlap = np.maximum(np.minimum(one[1], two[1])-np.maximum(one[0], two[0]), 0)
            if np.prod(overlap[plan]) < .70*np.prod(size[plan]):
                continue
            if abs(one.mean(axis=0)[axis]-two.mean(axis=0)[axis]) < .25*size[axis]:
                continue
            return {
                "compact_detail_solid_count": int(small),
                "shell_normal_axis_prototype": axis,
                "shell_bounding_dimensions_mm": _rounded([first[0], second[0]]),
                "shell_plan_overlap_ratio": round(float(np.prod(overlap[plan])/np.prod(size[plan])), 5),
                "shell_fill_ratios": [round(float(v/np.prod(d)), 5) for d, _, v in (first, second)],
            }
    return None


def _prototype(model, definition):
    data = model.features[definition]
    if not data.get("is_valid"):
        return None
    short, mid, long = sorted(data["dims"])
    bodies = data["solid_count"]
    if not 0 < short <= long:
        return None
    def measure(label, evaluator, with_data=False):
        shape = model.shape(definition)
        box = bounds(shape)
        result = evaluator(shape, box, data) if with_data else evaluator(shape, box)
        return (label, result) if result else None

    if bodies >= 20 and short <= .12*mid and mid >= 45:
        return measure("mesh", _mesh)
    if bodies >= 12 and 45 <= short and long <= 350 and long/short <= 2.75:
        return measure("alarm", _alarm_package)
    if bodies != 1:
        return None
    if (12 <= short <= 65 and 30 <= mid <= 95 and 40 <= long <= 140
            and long < 2.35*short + 20 and 25 <= data["faces"] <= 130):
        ring = measure("ring", _ring, True)
        if ring:
            return ring
    if (7 <= short <= 35 and long <= 90 and long >= 2.5*mid
            and len(data["cylinders"]) >= 6):
        return measure("clip", _clip, True)
    return None


def _neighbors(model, ident, box, ids, boxes, subtype):
    aabb_mm = np.linalg.norm(
        np.maximum(np.maximum(box[0]-boxes[:, 1], boxes[:, 0]-box[1]), 0),
        axis=1,
    )
    ordered = np.argsort(aabb_mm, kind="stable")
    minimum = next((float(aabb_mm[index]) for index in ordered if ids[index] != ident), None)
    nearby = [
        index for index in ordered if ids[index] != ident and aabb_mm[index] <= 4.0
    ]
    exact = []
    failed = []
    # Multi-body packages can hold hundreds of exact faces. Measure a physical
    # solid near each neighbor as a witness. Positive contact on that solid
    # applies to the package; a positive gap is local to the measured solid.
    source = model.local_shape(ident) if nearby else None
    probes = []
    if nearby and subtype in ("mesh", "alarm"):
        probes = [(solid, bounds(solid)) for solid in indexed(source, TopAbs_SOLID)]
    for index in nearby[:4]:
        try:
            if probes:
                target = boxes[index]
                def priority(probe):
                    limits = probe[1]
                    low = np.maximum(np.maximum(target[0]-limits[1], limits[0]-target[1]), 0)
                    return float(low@low)
                sample = min(probes, key=priority)[0]
            else:
                sample = source
            result = gap(sample, model.local_shape(ids[index]))
        except (AssertionError, RuntimeError, ValueError):
            failed.append(ids[index])
            continue
        exact.append({
            "id": ids[index],
            "box_distance_mm": round(float(aabb_mm[index]), 5),
            "gap_mm": round(result["distance_mm"], 5),
            "points": _rounded([result["point_1_mm"], result["point_2_mm"]]),
            "inner": bool(result["inner_solution"]),
            "geometry_scope": "one_solid" if probes else "whole_occurrence",
        })
    return sorted(exact, key=lambda item: item["gap_mm"]), len(nearby), failed, minimum


def measure_accessories(model, base):
    """Return candidate evidence with neutral occurrence IDs and no source names."""
    model.origin = np.asarray(base["basis"]["origin"])
    model.axes = np.asarray(base["basis"]["axes"])
    by_id = {model.ids[row["path"]]: row["definition"] for row in model.rows}
    parts = base["parts"]
    frame = np.asarray(base["frame_box"])
    frame_span = np.ptp(frame, axis=0)
    prototypes = {}
    for definition in set(by_id.values()):
        prototypes[definition] = _prototype(model, definition)
    candidates = {ident for ident, definition in by_id.items() if prototypes[definition]}
    ids = sorted(ident for ident, item in parts.items()
                 if "box" in item and item["valid"] and item["solids"] == 1
                 and item["volume"] >= 3000 and ident not in candidates)
    boxes = np.asarray([parts[ident]["box"] for ident in ids]).reshape(-1, 2, 3)
    collected = []
    for ident in sorted(candidates):
        subtype, metrics = prototypes[by_id[ident]]
        box = np.asarray(parts[ident]["box"])
        normalized_center = (box.mean(axis=0)-frame[0])/frame_span
        distance_to_sides = abs(np.minimum(normalized_center, 1-normalized_center))
        if subtype == "alarm" and (distance_to_sides[0] > .10 or distance_to_sides[1] > .18):
            continue
        if subtype == "mesh" and min(distance_to_sides[:2]) > .15:
            continue
        if subtype in ("ring", "clip") and min(distance_to_sides[:2]) > .22:
            continue
        exact, count, failed, minimum = _neighbors(model, ident, box, ids, boxes, subtype)
        receiver = next((item for item in exact if item["gap_mm"] <= .10), None)
        if receiver is None:
            receiver = next((item for item in exact if item["gap_mm"] <= 4.0), None)
        collected.append({
            "a": ident,
            "b": receiver["id"] if receiver else ident,
            "kind": "accessory_geometry",
            "subtype": subtype,
            "metrics": metrics,
            "box_mm": _rounded(box),
            "normalized_center": _rounded(normalized_center),
            "neighbor_screen_cap": 4,
            "neighbor_screen_radius_mm": 4.0,
            "screened_neighbor_count": count,
            "nearest_screened_population_box_gap_mm": round(minimum, 5) if minimum is not None else None,
            "neighbor_distance_failures": failed,
            "neighbor_witnesses": exact,
        })
    return collected
