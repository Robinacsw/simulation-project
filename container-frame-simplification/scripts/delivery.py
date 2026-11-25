"""Match approved geometry to the current CAD copy and check every exported body."""

from collections import Counter, defaultdict
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.spatial import cKDTree

from common import checked, load, save, sha, write_csv


HERE = Path(__file__).resolve().parent
MATCH_MM = .01
MATCH_V_REL = .00005


def center(row):
    xyz = row["center"]
    return np.asarray(xyz) if xyz is not None else np.asarray(row["tight_box"]).reshape(2,3).mean(axis=0)


def delta(a, b):
    return (float(np.max(np.abs(np.asarray(a["tight_box"])-b["tight_box"]))),
            float(np.linalg.norm(center(a)-center(b))),abs(a["volume"]-b["volume"]))


def matches(a, b, bridge=False):
    if (a["center"] is None) != (b["center"] is None):
        return False
    same_faces = a["faces"] == b["faces"]
    if bridge and "area" in a and "area" in b:
        # A periodic seam may split a face in STEP and leave it whole in X_T.
        # Require a fresh surface-area witness in addition to the shape/pose signature.
        if abs(a["faces"]-b["faces"])>1 or abs(a["area"]-b["area"])>max(.005,MATCH_V_REL*a["area"]):
            return False
    elif not same_faces:
        return False
    box, midpoint, volume = delta(a,b)
    return max(box,midpoint)<=MATCH_MM and volume<=max(.005,MATCH_V_REL*abs(a["volume"]))


class Inventory:
    def __init__(self, rows):
        self.rows = rows
        self.tree = cKDTree(np.asarray([center(r) for r in rows]).reshape((-1,3)))

    def options(self, wanted, bridge=False):
        return sorted(i+1 for i in self.tree.query_ball_point(center(wanted),MATCH_MM) if matches(wanted,self.rows[i],bridge))


def pair_subset(expected, observed):
    future = Inventory(observed)
    options = [future.options(row) for row in expected]
    edge_rows = np.repeat(np.arange(len(expected)),[len(row) for row in options])
    edge_cols = [i-1 for row in options for i in row]
    adjacency = csr_matrix((np.ones(len(edge_cols),dtype=bool),(edge_rows,edge_cols)),shape=(len(expected),len(observed)))
    pairing = maximum_bipartite_matching(adjacency,perm_type="column")
    missing = [i+1 for i,k in enumerate(pairing) if k < 0]
    if missing:
        raise ValueError(f"Cannot pair {len(missing)} expected bodies, first unmatched {missing[:12]}")
    return {i+1:int(k)+1 for i,k in enumerate(pairing)}


def pair_inventory(expected, observed):
    if len(expected) != len(observed):
        raise ValueError(f"Inventory count mismatch {len(expected)} / {len(observed)}")
    return pair_subset(expected,observed)


def select_occurrences(inventory, targets, protected=()):
    tree = Inventory(inventory)
    excluded = set(protected)
    assigned, skipped, claims = {}, {}, defaultdict(list)
    for ident,group in targets.items():
        choices = [tree.options(body,bridge=True) for body in group]
        if not choices or any(len(c)!=1 for c in choices):
            skipped[ident] = "当前几何无唯一实体；候选序号 " + str(choices)
            continue
        result = [c[0] for c in choices]
        if excluded.intersection(result) or len(set(result)) != len(group):
            skipped[ident] = "命中保护范围或一件多体对应发生重叠"
            continue
        assigned[ident] = result
        for body in result:
            claims[body].append(ident)
    for ident in {ident for ids in claims.values() if len(ids)>1 for ident in ids}:
        assigned.pop(ident)
        skipped[ident] = "多个匿名实例竞争同一实体"
    return assigned,skipped


def cad_job(spaceclaim, folder, source, mode, plan=None, target=None):
    folder.mkdir()
    request = folder/"request.json"
    save(request,{"source":str(source),"mode":mode,"target":str(target) if target else None,"plan":plan})
    env = dict(os.environ,FRAME_SIMPLIFICATION_REQUEST=str(request))
    argv = [str(spaceclaim),f"/RunScript={HERE/'spaceclaim_geometry.py'}","/ScriptAPI=252","/Headless=True",
            "/ExitAfterScript=True","/CommandLineAddInsOnly=true"]
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    with (folder/"stdout.log").open("w",encoding="utf-8") as out,(folder/"stderr.log").open("w",encoding="utf-8") as err:
        result = subprocess.run(argv,cwd=folder,env=env,stdout=out,stderr=err,creationflags=flags,timeout=1200)
    log = folder/"spaceclaim.log"
    if result.returncode or not log.exists() or not log.read_text(encoding="utf-8").rstrip().endswith("COMPLETE"):
        raise RuntimeError(f"CAD geometry {mode} failed, inspect {folder}")
    receipt = load(folder/"inventory.json")
    if receipt["count"] != len(receipt["rows"]):
        raise ValueError("Invalid CAD geometry inventory")
    print(f"CAD {mode}: {receipt['count']} bodies, {receipt['solids']} solid occurrences",flush=True)
    return receipt


def prepare_targets(analysis, approval):
    from intake import Model

    _,_,step,raw = analysis
    if raw["mode"] != "G+C":
        raise ValueError("Deliver requires connection-aware recognition")
    rules = set(approval["rules"])
    ids = set(approval.get("instances",[]))
    if not rules or not approval["decision"].strip():
        raise ValueError("Provide the engineering cleanup approval")
    if ids-raw["results"].keys():
        raise ValueError("Approval references unknown anonymous instances")
    chosen = {i:row for i,row in raw["results"].items()
              if row["rule"] in rules and (not ids or i in ids) and row["action"]=="清理几何"
              and (row["evidence"] or i in ids)}
    if ids and ids!=chosen.keys():
        raise ValueError("Some explicitly selected instances still lack supported clearance")
    model = Model(step)
    if model.paths.keys()!=raw["results"].keys():
        raise ValueError("STEP occurrence IDs cannot be reproduced")
    private = load(analysis[0]/"trace_private.json")
    if any(private[i]["occurrence"]!=model.paths[i] for i in chosen):
        raise ValueError("STEP occurrence order drifted since anonymization")
    targets = {i:model.signatures(i,require_complete=True) for i in chosen}
    # All structural keep witnesses with a verified position also protect the
    # actual CAD occurrence, even when an approved accessory appears coincident.
    retained = [signature for i,row in raw["results"].items() if row["action"]=="保留几何"
                for signature in model.signatures(i)]
    return chosen, targets, retained


def deliver(analysis_folder, geometry, approval_path, protection_path, spaceclaim, folder, name):
    from analysis_bundle import read_analysis

    folder = Path(folder).resolve()
    geometry = Path(geometry).resolve(strict=True)
    spaceclaim = Path(spaceclaim).resolve(strict=True)
    if Path(name).name!=name or Path(name).suffix.lower()!=".x_t" or geometry.suffix.lower() not in {".x_t",".step",".stp"}:
        raise ValueError("Use X_T/STEP input and a new simple .x_t output name")
    approval = load(approval_path)
    digest = approval["source_sha256"]
    checked(geometry,digest)
    external = load(protection_path) if protection_path else {"active":[],"other":[]}
    analysis = read_analysis(analysis_folder)
    folder.mkdir(parents=True)
    source = folder/("frozen_input"+geometry.suffix.lower())
    shutil.copyfile(geometry,source)
    checked(source,digest)
    baseline = cad_job(spaceclaim,folder/"read_before",source,"read")
    before = baseline["rows"]
    chosen,targets,keeps = prepare_targets(analysis,approval)
    tree = Inventory(before)
    active = external["active"]
    mandatory = active+external["other"]
    verified_protected = set(pair_subset(mandatory,before).values()) if mandatory else set()
    protected = {i for witness in mandatory for i in tree.options(witness)}
    protected.update(i for witness in keeps for i in tree.options(witness,bridge=True))
    selected, skipped = select_occurrences(before,targets,protected)
    removed = {i for ids in selected.values() for i in ids}
    staged = folder/"candidate.x_t"
    if removed or source.suffix.lower()!=".x_t":
        plan = {"source_bodies":baseline["count"],"source_solids":baseline["solids"],
                "removed_indices":sorted(removed),"protected_body_indices":sorted(protected),
                "before_samples":{str(i):before[i-1] for i in removed|protected}}
        checked(source,digest)
        written = cad_job(spaceclaim,folder/"write",source,"write",plan=plan,target=staged)
        pair_inventory(before,written["rows"])
    else:
        shutil.copyfile(source,staged)
    readback = cad_job(spaceclaim,folder/"read_after",staged,"read")
    expected = [r for i,r in enumerate(before,1) if i not in removed]
    compared = pair_inventory(expected,readback["rows"])
    if readback["solids"] != baseline["solids"]-len(removed) or readback["count"]-readback["solids"] != baseline["count"]-baseline["solids"]:
        raise ValueError("Unexpected solid or sheet count after independent readback")
    if mandatory:
        pair_subset(mandatory,readback["rows"])
    checked(geometry,digest)
    with staged.open("rb") as stream:
        header = stream.read(1024)
    if b"PARASOLID" not in header or b"FORMAT=text" not in header:
        raise ValueError("Output is not Parasolid text")
    final = folder/name
    if final != staged:
        staged.rename(final)
    fields = ["匿名实例","角色","判据","主体序号","导出动作","核对依据"]
    rows = [dict(zip(fields,(i,row["role"],row["rule"],",".join(map(str,selected.get(i,[]))),
             "已清理" if i in selected else "待复核",skipped.get(i,row["basis"])))) for i,row in chosen.items()]
    write_csv(folder/"清理前后对账.csv",fields,rows)
    deltas = [delta(expected[i-1],readback["rows"][j-1]) for i,j in compared.items()]
    receipt = {"file":str(final),"sha256":sha(final),"geometry_source":str(geometry),"source_sha256":digest,
               "analysis_sha256":sha(analysis[0]/"sealed_prediction.json"),"approval":approval["decision"],
               "before_solids":baseline["solids"],"after_solids":readback["solids"],"remaining_sheets":readback["count"]-readback["solids"],
               "removed":len(removed),"removed_by_rule":dict(Counter(chosen[i]["rule"] for i in selected)),"reconciled":len(expected),
               "cross_format_face_partitions_verified":sum(sig["faces"]!=before[index-1]["faces"]
                    for ident, indices in selected.items() for sig,index in zip(targets[ident],indices)),
               "protected_verified":len(verified_protected),"active_connections_verified":len(active),
               "protected_geometry_matches":len(protected),"unresolved_approved_candidates":skipped,
               "geometry_deltas":list(map(max,zip(*deltas))) if deltas else [],"mass_objects_added":0,"solver_started":False}
    save(folder/"X_T几何验收.json",receipt)
    print(receipt,flush=True)
