"""Compare a completed skill run against independently frozen project evidence."""

import argparse
from collections import Counter
from pathlib import Path

from PIL import Image
from analysis_bundle import read_analysis
from common import checked, load, save
from delivery import pair_inventory, pair_subset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("reference", "review-analysis", "approved-analysis", "review-images", "delivery", "out"):
        parser.add_argument(f"--{arg}",type=Path,required=True)
    args = parser.parse_args()
    ref = load(args.reference)
    _,early_manifest,_,early = read_analysis(args.review_analysis)
    _,late_manifest,_,late = read_analysis(args.approved_analysis)
    assert early_manifest["source_sha256"]==late_manifest["source_sha256"]==ref["step_sha256"]
    assert len(early["results"])==len(late["results"])==ref["expected_leaves"]
    for raw in (early,late):
        iso = raw["isolation"]
        assert iso["opened"]==["prompt.md","input.json"] and iso["knowledge_base_entries_visible"]==0
    assert all(late["results"][i]["action"]=="保留几何" for i,r in early["results"].items() if r["action"]=="保留几何")
    promotions = Counter(r["rule"] for i,r in late["results"].items()
                         if r["action"]=="清理几何" and early["results"][i]["action"]=="暂留待判")
    assert promotions==ref["approved_replay_promotions"]
    recognized_connections = Counter(r["rule"] for r in late["results"].values() if r["action"]=="保留几何")
    assert all(recognized_connections[k]==n for k,n in ref["expected_protection_rules"].items())
    review = load(args.review_images/"截图索引.json")
    pictured = 0
    for group in review:
        if group["image"]:
            with Image.open(args.review_images/group["image"]) as image:
                image.verify()
            pictured += group["count"]
    assert pictured>0
    assert {g["rule"]:g["count"] for g in review if g["image"]}==ref["expected_review_groups"]
    rec = load(args.delivery/"X_T几何验收.json")
    assert rec["removed_by_rule"]==ref["expected_roles"]
    assert rec["protected_verified"]==ref["protected_signatures"] and rec["active_connections_verified"]==ref["active_signatures"]
    assert not rec["unresolved_approved_candidates"] and rec["mass_objects_added"]==0 and not rec["solver_started"]
    checked(rec["file"],rec["sha256"])
    checked(ref["baseline"],rec["source_sha256"])
    before = load(args.delivery/"read_before/inventory.json")["rows"]
    after = load(args.delivery/"read_after/inventory.json")["rows"]
    reference = load(ref["reference_inventory"])["rows"]
    checked(ref["historical_x_t"],ref["historical_x_t_sha256"])
    pair_inventory(reference,after)
    surviving_baseline_ordinals = set(pair_subset(after,before).values())
    removed = [r for i,r in enumerate(before,1) if i not in surviving_baseline_ordinals]
    pair_inventory(ref["expected_removed_signatures"],removed)
    assert len(removed)==rec["removed"]
    report = {"assembly_instances":len(late["results"]),"rules":late["rules"],"review_candidates_pictured":pictured,
              "solids_before":rec["before_solids"],"solids_after":rec["after_solids"],"sheets":rec["remaining_sheets"],
              "removals_by_role":rec["removed_by_rule"],"historical_reference_bodies_matched":len(after),
              "replay_promotions":dict(promotions),
              "anonymous_protection_counts":{k:recognized_connections[k] for k in ref["expected_protection_rules"]},
              "exact_delta_bodies_matched":len(removed),"active_connections_verified":rec["active_connections_verified"],
              "protected_occurrences_verified":rec["protected_verified"],"mass_objects_added":rec["mass_objects_added"],
              "solver_started":rec["solver_started"],"source_step_sha256":ref["step_sha256"],"export_receipt":rec["file"],
              "file_sha256":rec["sha256"],"picked_candidate_truth_used_before_seal":False}
    save(args.out,report)
    print(report)


if __name__ == "__main__":
    main()
