"""Measure side-door pivot and end-cladding relations from anonymous geometry."""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations
from math import dist

import numpy as np

from measure_geometry import gap


def _cylinder_reader(model):
    @lru_cache(maxsize=None)
    def projected(ident):
        # Model resolves the B-Rep occurrence by neutral ID. The assembly path
        # belongs to the adapter and supplies no decision feature.
        occurrence = model.by_path[model.paths[ident]]
        cylinders = model.features[occurrence["definition"]]["cylinders"]
        rotation, offset = model.pose(ident)
        local_rotation = model.axes @ rotation
        local_offset = model.axes @ (offset - model.origin)
        for_axis = []
        for cylinder in cylinders:
            start = local_rotation @ cylinder["a"] + local_offset
            end = local_rotation @ cylinder["b"] + local_offset
            length = np.linalg.norm(end - start)
            if length > .1 and abs(end[2] - start[2]) >= .985 * length:
                xy = (start[:2] + end[:2]) / 2
                displacement = np.linalg.norm(start[:2] - end[:2])
                if displacement <= .05:
                    for_axis.append((xy, sorted((start[2], end[2])), cylinder))
        return for_axis

    return projected


def _full_cylinders(cylinders):
    """Join angular patches sharing an axis, radius and axial interval."""
    angles = defaultdict(float)
    for xy, interval, cylinder in cylinders:
        if cylinder["angle"] >= 1.4:
            key = (tuple(np.round(xy, 3)), tuple(np.round(interval, 3)),
                   round(cylinder["r"], 3), cylinder["inner"])
            angles[key] += cylinder["angle"]
    return [key for key, angle in angles.items() if angle >= 6.25]


def _ear_pairs(cylinders, box):
    height = box[1][2] - box[0][2]
    ears = defaultdict(list)
    for axis, (low, high), radius, inner in cylinders:
        if inner and 3 <= radius <= 20 and .035 * height <= high - low <= .23 * height:
            ears[axis, radius].append((low, high))
    for (axis, radius), intervals in ears.items():
        intervals.sort()
        lower, upper = intervals[0], intervals[-1]
        if (upper[0] - lower[1] >= .34 * height
                and abs(lower[0] - box[0][2]) <= .2 * height
                and abs(upper[1] - box[1][2]) <= .2 * height):
            yield axis, radius, lower, upper


def _pivot_witness(axis, radius, lower, upper, neighbors):
    pins = []
    barrels = []
    for ident, sections in neighbors:
        for center, (low, high), r, inner in sections:
            offset = dist(axis, center)
            if offset > .08:
                continue
            if (not inner and .01 <= radius - r <= 1.2
                    and low <= lower[0] + .05 and high >= upper[1] - .05):
                pins.append((offset, ident, r, low, high))
            if (inner and radius - .3 <= r <= radius + 5
                    and lower[1] - .1 <= low < high <= upper[0] + .1
                    and high - low >= .27 * (upper[0] - lower[1])):
                barrels.append((offset, ident, low, high))
    pairs = [(pin, barrel) for pin in pins for barrel in barrels if pin[1] != barrel[1]]
    if not pairs:
        return None
    pin, barrel = min(pairs, key=lambda pair: (pair[0][0] + pair[1][0], pair[0][1], pair[1][1]))
    return {
        "axis_xy_mm": list(axis), "bore_radius_mm": radius,
        "ear_intervals_mm": [list(lower), list(upper)],
        "pin": pin[1], "pin_radius_mm": pin[2], "pin_interval_mm": list(pin[3:]),
        "barrel": barrel[1], "barrel_interval_mm": list(barrel[2:]),
        "pin_center_offset_mm": pin[0], "barrel_center_offset_mm": barrel[0],
    }


def _measure_side_pivots(parts, groups, frame, projected):
    span = np.ptp(frame, axis=0)

    @lru_cache(maxsize=None)
    def full(ident):
        return _full_cylinders(projected(ident))

    for ident, part in parts.items():
        box = np.asarray(part["box"])
        x, y, _ = (box.mean(axis=0) - frame[0]) / span
        small, middle, large = part["shape_mm"]
        if not (part["solids"] == 1 and .055 <= x <= .945 and min(y, 1-y) <= .05
                and 22 <= small <= 75 and 35 <= middle <= 115 and 45 <= large <= 175):
            continue
        pairs = list(_ear_pairs(full(ident), box))
        if not pairs:
            continue
        sibling_geometry = [(other, full(other)) for other in groups[part["parent"]]
                            if other != ident and parts[other]["shape_mm"][-1] <= 2.2 * large]
        matches = [_pivot_witness(*pair, sibling_geometry) for pair in pairs]
        if all(matches):
            yield {"a": ident, "b": matches[0]["pin"], "kind": "accessory_geometry",
                   "subtype": "side_hinge", "source": "position and B-Rep ears, shaft and hinge barrel",
                   "side_fraction": float(y), "pivot_witnesses": matches}


def _corrugation_profile(part, projected):
    depth, width, height = np.ptp(part["box"], axis=0)
    folds = np.asarray([xy for xy, (low, high), cylinder in projected
                        if .1 <= cylinder["angle"] < 1.4 and cylinder["r"] <= .33 * depth
                        and high - low >= .065 * height])
    if len(folds) < 6:
        return None
    turns = sorted(set(np.round(folds[:, 1], 1)))
    fold_depth = float(np.ptp(folds[:, 0]))
    if len(turns) >= 6 and turns[-1] - turns[0] >= .53 * width and fold_depth >= .33 * depth:
        return {"folds_along_width_mm": turns, "bend_axis_depth_range_mm": fold_depth,
                "thickness_proxy_mm": 2 * part["volume"] / part["area"]}
    return None


def _near_end_rails(ident, parts, frame, test_gap):
    panel = parts[ident]
    box = np.asarray(panel["box"])
    panel_width = box[1, 1] - box[0, 1]
    frame_width = frame[1][1] - frame[0][1]
    candidates = [[], []]
    for other, host in parts.items():
        if host["parent"] != panel["parent"] or host["solids"] != 1 or other == ident:
            continue
        lo, hi = np.asarray(host["box"])
        span = hi - lo
        if not (span[1] >= .60 * frame_width and min(span[0], span[2]) >= 20
                and 3 * max(span[0], span[2]) < span[1]):
            continue
        if min(box[1, 1], hi[1]) - max(box[0, 1], lo[1]) < .65 * panel_width:
            continue
        if max(box[0, 0] - hi[0], lo[0] - box[1, 0]) > 12:
            continue
        top = (lo[2] + hi[2]) / 2 > box[:, 2].mean()
        if abs((lo[2] - box[1, 2]) if top else (hi[2] - box[0, 2])) > 12:
            continue
        measured = test_gap(ident, other)
        if measured is not None and measured["distance_mm"] <= 12 and not measured["inner_solution"]:
            candidates[int(top)].append({"member": other, **measured})
    if all(candidates):
        return [min(choices, key=lambda row: (row["distance_mm"], row["member"])) for choices in candidates]
    return None


def _measure_end_cladding(model, parts, frame, projected):
    span = np.ptp(frame, axis=0)
    panels = {}
    for ident, part in parts.items():
        box = np.asarray(part["box"])
        depth, width, height = np.ptp(box, axis=0)
        x = (box[:, 0].mean() - frame[0][0]) / span[0]
        if not (part["solids"] == 1 and min(x, 1-x) <= .055 and height >= .56 * span[2]
                and width >= .19 * span[1] and 12 <= depth <= .13 * width
                and 0 < part["area"] and 2 * part["volume"] / part["area"] <= .35 * depth):
            continue
        profile = _corrugation_profile(part, projected(ident))
        if profile:
            panels[ident] = (box, profile)

    @lru_cache(maxsize=None)
    def test_gap(ident, other):
        try:
            return gap(model.local_shape(ident), model.local_shape(other))
        except (AssertionError, RuntimeError, ValueError):
            return None

    for a, b in combinations(panels, 2):
        if parts[a]["parent"] != parts[b]["parent"]:
            continue
        box_a, profile_a = panels[a]
        box_b, profile_b = panels[b]
        seam = max(box_a[0, 1] - box_b[1, 1], box_b[0, 1] - box_a[1, 1])
        coverage = max(box_a[1, 1], box_b[1, 1]) - min(box_a[0, 1], box_b[0, 1])
        if not (-.10 * min(np.ptp(box_a[:, 1]), np.ptp(box_b[:, 1])) <= seam <= 12
                and coverage >= .60 * span[1] and abs(box_a[:, 0].mean() - box_b[:, 0].mean()) <= 12
                and np.max(np.abs(box_a[:, 2] - box_b[:, 2])) <= .02 * span[2]):
            continue
        rail_a = _near_end_rails(a, parts, frame, test_gap)
        rail_b = _near_end_rails(b, parts, frame, test_gap)
        if rail_a and rail_b and [r["member"] for r in rail_a] == [r["member"] for r in rail_b]:
            for ident, other, profile, rails in ((a, b, profile_a, rail_a), (b, a, profile_b, rail_b)):
                yield {"a": ident, "b": other, "kind": "accessory_geometry", "subtype": "end_corrugation",
                       "source": "paired corrugated sheet and B-Rep gaps to end rails", **profile,
                       "seam_envelope_gap_mm": float(seam), "end_member_witnesses": rails}


def measure_side_geometry(model, base):
    """Use prepared occurrence poses and B-Rep data, anonymous hierarchy and box axes."""
    parts = {ident: part for ident, part in base["parts"].items() if "box" in part and part["valid"]}
    frame = np.asarray(base["frame_box"])
    groups = defaultdict(list)
    for ident, part in parts.items():
        groups[part["parent"]].append(ident)
    projected = _cylinder_reader(model)
    edges = list(_measure_side_pivots(parts, groups, frame, projected))
    edges.extend(_measure_end_cladding(model, parts, frame, projected))
    return edges
