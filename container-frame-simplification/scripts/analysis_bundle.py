"""Prepare only anonymous CAD evidence for the rules subprocess, then seal the result."""

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import shutil
import stat
import subprocess
import sys

from common import checked, load, save, sha, write_csv


HERE = Path(__file__).resolve().parent
PROMPT = HERE.parent / "assets/anonymous-task.md"
REVISION = "frame-EU-1.5"
STACK = ("shaft_stack", "axial_gate", "accessories", "roof", "roof_holder")
APPROVAL_FLAGS = {"roof_connector_approved", "local_roof_holder_approved", "bottom_pad_underskin_approved"}


def profile(path):
    source = load(path)
    if set(source) != {"cases", "scope"}:
        raise ValueError("Profile must contain just cases and scope")
    fields = {"name", "entrance", "support"}
    cases = source["cases"]
    if (not isinstance(cases,list) or not cases or
        any(set(c) != fields or any(not isinstance(v,str) or not v.strip() for v in c.values()) for c in cases)):
        raise ValueError("Supply each selected test case name, load entrance and support")
    scope = source["scope"]
    boundary = scope.get("fixture_boundary")
    if boundary not in {"discrete_castings_and_floor", "other", "unconfirmed"}:
        raise ValueError("Set the fixture boundary explicitly")
    if set(scope) - APPROVAL_FLAGS - {"fixture_boundary"} or any(type(scope[k]) is not bool for k in APPROVAL_FLAGS & scope.keys()):
        raise ValueError("Use only the documented anonymous scope flags")
    return cases, scope


def result_csv(folder, result, trace):
    fields = ["匿名实例", "原始标签", "角色", "建议动作", "判据", "证据", "待核条件"]
    rows = [dict(zip(fields, (ident, " / ".join(trace[ident]["labels"]), row["role"], row["action"],
             row["rule"], ",".join(row["evidence"]), row["falsifier"]))) for ident,row in result["results"].items()]
    write_csv(folder/"逐件识别与去留.csv",fields,rows)


def run_worker(folder, base, relations, mode, case_profile):
    from recognize import RULE_VERSION

    cases,scope = profile(case_profile)
    pieces = [HERE/"recognize.py"]
    if mode == "G+C":
        pieces.extend(HERE/"rules"/(name+".py") for name in STACK)
    bundle = b"\n".join(p.read_bytes() for p in pieces)
    worker = folder/"anonymous_worker"
    worker.mkdir()
    prompt = worker/"prompt.md"
    shutil.copyfile(PROMPT,prompt)
    packet = worker/"input.json"
    save(packet,{"mode":mode,"base":base,"rules":REVISION if mode=="G+C" else RULE_VERSION,
                 "authorization":"只读本任务包","cases":cases,"scope":scope,"relations":relations if mode=="G+C" else None,
                 "run":"anonymous-run"})
    result = subprocess.run([sys.executable,"-I","-S","-B","-c",(HERE/"worker_code.txt").read_text(encoding="utf-8")],
                            cwd=worker,input=bundle,capture_output=True,timeout=180)
    if result.returncode:
        raise RuntimeError(result.stderr.decode("utf-8",errors="replace"))
    report = json.loads(result.stdout)
    iso = report["isolation"]
    if (iso["opened"] != ["prompt.md","input.json"] or iso["denied_probes"] != ["open","os.scandir"]
        or iso["knowledge_base_entries_visible"] or iso["retrieval_tools"]):
        raise ValueError("Anonymous subprocess read boundary not verified")
    if (iso["input_sha256"],iso["rules_sha256"],report["prompt_sha256"]) != (sha(packet),sha256(bundle).hexdigest().upper(),sha(prompt)):
        raise ValueError("Anonymous run fingerprints did not match")
    save(folder/"sealed_prediction.json",report)
    for path in (packet,prompt,folder/"sealed_prediction.json"):
        path.chmod(stat.S_IREAD)
    print(f"{mode} {report['rules']}: {report['actions']}",flush=True)
    return report


def measure_connections(model, base):
    from accessory_geometry import measure_accessories
    from measure import Scanner
    from relations import Graph
    from roof_geometry import measure_roof_geometry
    from side_geometry import measure_side_geometry
    from support_geometry import measure_support_geometry

    graph = Graph(model,base)
    prior = graph.run()["edges"]
    joint = Scanner(graph).run(prior)
    print("Clamping and facing surfaces",Counter(e["kind"] for e in joint),flush=True)
    extra = []
    for collector in (measure_accessories,measure_side_geometry,measure_support_geometry,measure_roof_geometry):
        measured = collector(model,base)
        print(collector.__name__,Counter(e["subtype"] for e in measured),flush=True)
        extra.extend(measured)
    for n,edge in enumerate(sorted(extra,key=lambda e:(e["subtype"],e["a"])),1):
        edge.setdefault("b",edge["a"])
        edge["edge_id"] = f"R{n:05d}"
    return prior+joint+extra


def scan_assembly(step, configuration, case_profile, folder, mode, expected_leaves=None, datum=None):
    from intake import Model

    folder = Path(folder).resolve()
    folder.mkdir(parents=True)
    profile(case_profile)
    model = Model(step)
    if expected_leaves is not None and len(model.rows) != expected_leaves:
        raise ValueError(f"STEP leaf census expected {expected_leaves}, read {len(model.rows)}")
    model.prepare_definitions()
    base = model.make_base(load(datum) if datum else None)
    relations = measure_connections(model,base) if mode == "G+C" else []
    trace = model.private_trace()
    for name,value in [("base.json",base),("relations.json",relations),("trace_private.json",trace)]:
        save(folder/name,value)
    report = run_worker(folder,base,relations,mode,case_profile)
    result_csv(folder,report,trace)
    finish(folder,model.source,model.source_hash,configuration,report)


def finish(folder, source, fingerprint, configuration, report):
    files = ("base.json","relations.json","trace_private.json","sealed_prediction.json")
    save(folder/"manifest.json",{"source":str(source),"source_sha256":fingerprint,"configuration":configuration,
         "mode":report["mode"],"rules":report["rules"],"files":{p:sha(folder/p) for p in files},
         "decoded_occurrences":len(report["results"])})


def read_analysis(folder):
    folder = Path(folder).resolve(strict=True)
    manifest = load(folder/"manifest.json")
    required = {"base.json","relations.json","trace_private.json","sealed_prediction.json"}
    if set(manifest["files"]) != required:
        raise ValueError("Incomplete analysis manifest")
    for name,expected in manifest["files"].items():
        checked(folder/name,expected)
    source = checked(manifest["source"],manifest["source_sha256"])
    sealed = load(folder/"sealed_prediction.json")
    return folder, manifest, source, sealed


def replay(prior, case_profile, folder, mode=None):
    prior,manifest,source,_ = read_analysis(prior)
    mode = manifest["mode"] if mode is None else mode
    if manifest["mode"] == "G" and mode == "G+C":
        raise ValueError("Run G+C with CAD measurement to add connections")
    folder = Path(folder).resolve()
    folder.mkdir(parents=True)
    for name in ("base.json","relations.json","trace_private.json"):
        shutil.copyfile(prior/name,folder/name)
    report = run_worker(folder,load(folder/"base.json"),load(folder/"relations.json"),mode,case_profile)
    result_csv(folder,report,load(folder/"trace_private.json"))
    finish(folder,source,manifest["source_sha256"],manifest["configuration"],report)
