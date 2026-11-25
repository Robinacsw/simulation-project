"""Measure geometric hosts for small assemblies along the lower side rails."""

import numpy as np

from measure import Scanner, TOUCH_MM
from relations import Graph, intervals


FRAME_ROLES = {"BOTTOM_RAIL", "BOTTOM_LONG", "CROSS", "POST", "CAST", "DECK"}
HOST_SEARCH_MM = 4.0


class PatchScanner(Scanner):
    """Use the existing planar-face probe with a neutral shape provider."""

    def __init__(self, graph):
        self.g = graph


def lower_candidates(graph):
    pads, angles, ribs = [], [], []
    middle_y = np.mean(graph.frame, axis=0)[1]
    for ident in graph.ids:
        index = graph.look[ident]
        if graph.lo[index, 2] > graph.frame[0][2] + .08 * graph.span[2]:
            continue
        if abs(graph.center[index, 1] - middle_y) < .37 * graph.span[1]:
            continue
        x, y, z = graph.size[index]
        part = graph.parts[ident]
        occupancy = part["volume"] / max(1.0, x * y * z)
        if 4 <= z <= .25 * min(x, y) and 30 <= min(x, y) <= max(x, y) <= 450:
            if part["faces"] == 6 and occupancy >= .99:
                pads.append(ident)
        if 55 <= x <= 450 and .30 * x <= min(y, z) and .60 <= y / max(1.0, z) <= 1.65:
            if .12 <= occupancy <= .48:
                angles.append(ident)
        if 2 <= x <= .24 * min(y, z) and 25 <= min(y, z) <= max(y, z) <= 240:
            if .40 <= occupancy <= .74 and .48 <= y / z <= 2.1:
                ribs.append(ident)
    return pads, angles, ribs


def search_near(graph, ident, candidates, distance):
    index = graph.look[ident]
    separation = intervals(graph.lo[candidates], graph.hi[candidates], graph.lo[index], graph.hi[index])
    return [graph.ids[j] for j in candidates[np.max(separation, axis=1) <= distance] if j != index]


def single_host(graph, scanner, ident, frame_indices):
    touching, nearby = [], []
    for other in search_near(graph, ident, frame_indices, HOST_SEARCH_MM):
        witness = graph.nearest(ident, other)
        if witness is None:
            return None
        if witness["gap_mm"] <= TOUCH_MM:
            touching.append((other, witness))
        elif witness["gap_mm"] <= HOST_SEARCH_MM:
            nearby.append({"instance": other, "role_shape": graph.types[other], "nearest": witness,
                           "patch": scanner.patch(ident, other, tolerance=HOST_SEARCH_MM)})
    if len(touching) != 1:
        return None
    host, witness = touching[0]
    if graph.types[host] == "BOTTOM_RAIL" and not witness["inner"]:
        return host, witness, nearby
    return None


def plane_axes(graph, ident, patch):
    number = patch["faces"][0]
    return next((np.abs(normal) for face, _, _, _, normal in graph.planes(ident) if face == number), None)


def record(graph, ident, backing, host, subtype, group, patch, witness, neighbors):
    index = graph.look[ident]
    unresolved_facing = any(other["patch"] is not None for other in neighbors)
    return {"a": ident, "b": backing, "kind": "accessory_geometry", "subtype": subtype,
            "host": host, "group": group, "patch": patch, "nearest": witness,
            "center_mm": graph.center[index].round(5).tolist(),
            "box_mm": graph.parts[ident]["box"], "structural_neighbor_search_mm": HOST_SEARCH_MM,
            "additional_frame_neighbors": neighbors, "weld_interface_review": unresolved_facing,
            "source": "B-Rep nearest points and opposed trimmed-plane overlap"}


def measure_support_geometry(model, base):
    graph = Graph(model, base)
    scanner = PatchScanner(graph)
    pads, angles, ribs = lower_candidates(graph)
    frame_indices = np.array([graph.look[i] for i in graph.ids if graph.types[i] in FRAME_ROLES], dtype=int)
    rail_indices = np.array([i for i in frame_indices if graph.types[graph.ids[i]] == "BOTTOM_RAIL"], dtype=int)
    rib_indices = np.array([graph.look[i] for i in ribs], dtype=int)
    middle_y = np.mean(graph.frame, axis=0)[1]
    edges = []

    for ident in pads:
        i = graph.look[ident]
        rail_options = search_near(graph, ident, rail_indices, TOUCH_MM)
        if not any(abs(graph.hi[i, 2] - graph.lo[graph.look[rail], 2]) <= TOUCH_MM for rail in rail_options):
            continue
        found = single_host(graph, scanner, ident, frame_indices)
        if found is None:
            continue
        host, witness, neighbors = found
        i, j = graph.look[ident], graph.look[host]
        if abs(graph.hi[i, 2] - graph.lo[j, 2]) > TOUCH_MM:
            continue
        patch = scanner.patch(ident, host)
        if patch and patch["projected_area_mm2"] >= .45 * graph.size[i, 0] * graph.size[i, 1]:
            edges.append(record(graph, ident, host, host, "bottom_pad", [ident], patch, witness, neighbors))

    for ident in angles:
        rib_options = search_near(graph, ident, rib_indices, TOUCH_MM)
        if len(rib_options) < 2:
            continue
        found = single_host(graph, scanner, ident, frame_indices)
        if found is None:
            continue
        host, witness, neighbors = found
        i, j = graph.look[ident], graph.look[host]
        if abs(graph.center[i, 1] - middle_y) <= abs(graph.center[j, 1] - middle_y):
            continue
        if abs(graph.lo[i, 2] - graph.lo[j, 2]) > .22 * graph.size[i, 2]:
            continue
        patch = scanner.patch(ident, host)
        if not patch or patch["projected_area_mm2"] < .35 * graph.size[i, 0] * graph.size[i, 2]:
            continue
        if plane_axes(graph, ident, patch)[1] < .999:
            continue
        # The other broad limb is parallel to the bottom plane.
        if not any(abs(normal[2]) >= .999 and area >= .35 * graph.size[i, 0] * graph.size[i, 1]
                   for _, _, area, _, normal in graph.planes(ident)):
            continue

        rib_witnesses = []
        for rib in rib_options:
            k = graph.look[rib]
            offset = graph.center[k, 0] - graph.center[i, 0]
            if not .28 <= abs(offset) / graph.size[i, 0] <= .60:
                continue
            if not any(abs(normal[0]) < .05 and min(abs(normal[1]), abs(normal[2])) > .35
                       for _, _, _, _, normal in graph.planes(rib)):
                continue
            nearby = graph.nearest(rib, ident)
            attachment = scanner.patch(rib, ident)
            if nearby and not nearby["inner"] and nearby["gap_mm"] <= TOUCH_MM and attachment:
                if attachment["projected_area_mm2"] >= .25 * graph.size[k, 0] * min(graph.size[k, 1:]):
                    rib_witnesses.append((offset, rib, attachment, nearby))

        if len(rib_witnesses) != 2 or rib_witnesses[0][0] * rib_witnesses[1][0] >= 0:
            continue
        other_hosts = []
        for _, rib, _, _ in rib_witnesses:
            for other in search_near(graph, rib, frame_indices, HOST_SEARCH_MM):
                nearby = graph.nearest(rib, other)
                if nearby is None or (other != host and nearby["gap_mm"] <= HOST_SEARCH_MM):
                    other_hosts.append(other)
        if other_hosts:
            continue
        group = [ident, *sorted(rib for _, rib, _, _ in rib_witnesses)]
        edges.append(record(graph, ident, host, host, "fixture_angle", group, patch, witness, neighbors))
        for _, rib, attachment, nearby in rib_witnesses:
            edges.append(record(graph, rib, ident, host, "fixture_gusset", group, attachment, nearby, []))

    return edges
